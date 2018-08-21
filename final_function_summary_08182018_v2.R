#!/usr/bin/env Rscript

##packages :
rm(list=ls())
library(plyr)
library(dplyr)
library('RSQLite')
library(grDevices)
library(splitstackshape)
library(R.utils)


input_args = commandArgs(trailingOnly=TRUE)

###Defining functions and arguments 

## Args 
args1 = input_args[1]

args4 = paste0(input_args[2],'/',input_args[4])

input_files_arg <- input_args[5]
#input_files_arg <- 'as_tmp.txt,as_tmp1.txt,as_tmp2.txt'
input_files_arg <- unlist(strsplit(input_files_arg,','))
input_files_arg5 <- paste0(input_args[2],'/',input_files_arg[1]) # this is used as the input file argument
input_files_arg5_ellipsis <- setdiff(input_files_arg,input_files_arg[1]) # this is the ... argument
input_files_arg5_ellipsis <- paste0(input_args[2],'/',input_files_arg5_ellipsis)

args6 = paste0(input_args[2],'/',input_args[6])
args7 = paste0(input_args[2],'/',input_args[7])


##Functions 

## Loads file with gene names to plot
# Input : Text file with gene names as rows
# Output: Table of 1 column with gene names

load_gene_list = function(fname){
  final_list = read.table(fname, stringsAsFactors = F)
  return(final_list)
}


## Makes SQL db called Test.sqlite 
# Input : large text file, exon, bed file, cosmic file 
# Output: table sin db with three tables 

load_largefile_to_sql = function(largefile,dblocation_and_name,new_load=T){
  db = dbConnect(SQLite(), dbname = dblocation_and_name)           
  #making the first table (tmp) in Test.sqlite from the text file
  if (new_load){
    dbSendQuery(db,'drop table if exists tmp')
    dbWriteTable(conn = db, name = 'tmp'                        
                 , value = largefile
                 , sep="\t", row.names = FALSE, header = FALSE, overwrite=T)
  } else {
    # create temporary table to hold data
    dbSendQuery(db,'drop table if exists tmp_new_data')
    dbWriteTable(conn = db, name = 'tmp_new_data'                        
                 , value = largefile
                 , sep="\t", row.names = FALSE, header = FALSE, overwrite=T)
    
    #insert the data from temporary table into the main "tmp" table
    dbSendQuery(db,'insert into tmp select * from tmp_new_data')
    
    # drop temporary table
    dbSendQuery(db,'drop table if exists tmp_new_data')
  }
  #diconnect 
  dbDisconnect(db)
  
}

load_file_to_sql = function(bedfile,cosmicfile,dblocation_and_name){
  db = dbConnect(SQLite(), dbname = dblocation_and_name)           
  
  #making the second table (tmp_bed) in Test.sqlite from the text.exon.bed file 
  dbSendQuery(db,'drop table if exists tmp_bed')            
  dbWriteTable(conn = db, name = 'tmp_bed'
               , value = bedfile
               , sep="\t", row.names = FALSE, header = FALSE, overwrite=T)
  
  #getting cosmic 
  dbSendQuery(db,'drop table if exists cosmic')            
  dbWriteTable(conn = db, name = 'cosmic'
               , value = cosmicfile
               , sep="\t", overwrite=T) 
  
  #diconnect 
  dbDisconnect(db)
}


###Querying the sql dbs


##Largetextfile 
#input: gene name (loop throughtthe genelist)
#output:text_file
small_data_1 = function(x, dblocation_and_name){                  
  db = dbConnect(SQLite(), dbname = dblocation_and_name) 
  querytext = paste0("select * from tmp where V1 = '", x ,"'")
  new_data_1 = dbGetQuery(db,querytext)
  colnames(new_data_1) = c('gene', 'chr','pos', 'OrigDepth', 'Depth_30', 'RefBase', "FwdSig", 'RevSig')
  new_data = ddply(new_data_1, c('gene','chr','pos','RefBase'), summarize,
                   Depth_30 = mean(as.numeric(Depth_30)), 
                   OrigDepth = mean(as.numeric(OrigDepth)), FwdSig = mean(as.numeric(FwdSig)),
                   RevSig = mean(as.numeric(RevSig)))
  dbDisconnect(db)
  return(new_data)
}


##Exonbedfile 
#input: gene name (loop throughtthe genelist)
#output:exon_bed_file
small_data_2 = function(x, dblocation_and_name){                               
  db_2 = dbConnect(SQLite(), dbname = dblocation_and_name) 
  querytext_2 = paste0("select * from tmp_bed where V4 = '", x ,"'")
  new_data_2 = dbGetQuery(db_2,querytext_2)
  dbDisconnect(db_2)
  return(new_data_2)
}


##calculating depth 
calculate_depth = function(Depth_file){
  gene_depth_files = data.frame() 
  Depth_file = as.data.frame(filter(Depth_file, Depth_30 < 100))
  if (nrow(Depth_file)>0){
    Depth_file = subset(Depth_file, select=-c(OrigDepth, end_pos, FwdSig, RevSig,transcript_position))
    names(Depth_file)[names(Depth_file)=="V5"] <- "Exons"
    gene_depth_files = rbind(gene_depth_files, Depth_file)
    regions = data.frame(seqToIntervals(gene_depth_files$trans_pos))
    gene_depth_files_region_1 = merge.data.frame(regions, gene_depth_files,
                                                 by.x = c("from"), 
                                                 by.y = c("trans_pos"), all.x = T)
    gene_depth_files_region_1$Region <- paste(gene_depth_files_region_1$from,"-", 
                                              gene_depth_files_region_1$to)
    gene_depth_files_region_1 = subset(gene_depth_files_region_1, select=-c(from, to))
    
    return(gene_depth_files_region_1)
  } else {
    return(paste0(unique(Depth_file$gene),' has good coverage'))
  }
}

## calculating directional noise 

calculate_directional_noise = function(temp_2){
  
  temp_2$dir_noise = abs(temp_2$FwdSig + temp_2$RevSig)
  temp_2 = filter(temp_2, dir_noise > 0.05)
  temp_2$dir_noise_norm = (abs(temp_2$FwdSig - temp_2$RevSig)/ abs(temp_2$FwdSig + temp_2$RevSig))
  temp_2$dir_noise_norm = sub(NaN, 0, temp_2$dir_noise_norm)
  temp_2 = filter(temp_2, dir_noise_norm > 0.5)
  names(temp_2)[names(temp_2)=="V5"] <- "Exons"
  temp_2 = filter(temp_2, Depth_30 >100)
  colnames(temp_2) = c("Position", "Exon", "Gene", "Chr", "End", "OrigDepth", "Depth30", "FwdSig", "RevSig", "TransPos", "Dir_noise", "Score")  
  temp_2 = temp_2[,c("Gene", "Chr","Position",  "Exon", "Score")]
  return(temp_2)
  
}  
##added today 
#Adding region 
#  adding_regions = function(genedirectionalnoise_template){
#  regions_2 = data.frame(seqToIntervals(genedirectionalnoise_template$trans_pos))
#  gene_directional_noise_merge_1 = merge.data.frame(regions_2, genedirectionalnoise_template, 
#                                                    by.x = c("from"), by.y = c("trans_pos"), all.x = T)
#  gene_directional_noise_merge_1$Region <- paste(gene_directional_noise_merge_1$from,"-", gene_directional_noise_merge_1$to)
#  gene_directional_noise_merge_1 = subset(gene_directional_noise_merge_1, 
#                                          select=-c(from, to, end_pos,OrigDepth, Depth_30, 
#                                                    FwdSig, RevSig,transcript_position, dir_noise))
#  return(gene_directional_noise_merge_1)
#}

total_bases_with_direct_noise = function(gene_directional_noise_merge_template, n){
  gene_directional_noise_merge_template$dir_noise_norm <- as.numeric(as.character(gene_directional_noise_merge_template$dir_noise_norm))
  gene_direct_noise_2 = nrow(filter(gene_directional_noise_merge_template, dir_noise_norm >= n))
  
  if (length(gene_direct_noise_2) > 0){
    return( gene_direct_noise_2)
  } else {
    return(0)
  }
  
}


##cosmic
#input: gene name (loop throughtthe genelist)
#output:cosmic_file3
small_data_3 = function(x){                  
  db_3 = dbConnect(SQLite(), dbname = 'Test.sqlite') 
  querytext_3 = paste0('select "Gene.name", "Mutation.Genome.Position", "SNP" 
                       from cosmic where "Gene.name" 
                       like "', x ,'%"  ')
  new_data3 = dbGetQuery(db_3,querytext_3)
  colnames(new_data3) = c("Gene.name", "Mutation.Genome.Position", "SNP")
  return(new_data3)
}


##

opening_exon_file = function(ebf_template){
  df2 = data.frame()
  for (i in 1:nrow(ebf_template)){
    
    col1 = c(ebf_template[i,'V2']:ebf_template[i,'V3'])
    df=as.data.frame(
      cbind(trans_pos = col1,
            gene = rep(ebf_template[i,'V4'],length(col1)),
            chr = rep(ebf_template[i,'V1'],length(col1)),
            end_pos = rep(ebf_template[i,'V3'],length(col1)),
            exon = rep(ebf_template[i,'V5'], length(col1))),
      stringsAsFactors = FALSE)
    colnames(df) = c('trans_pos', 'gene','chr', 'end_pos', 'exon')
    df2 = rbind(df2, df)
    
  }
  return(df2)
}



## combinign largetextfile & exon_bed files 
#input : text_file, exon bed file after opening it) 
#output : text_exon_bed

make_text_exon_bed = function(template_file, input_file){
  get_text_exon_bed = merge(template_file[ ,c('trans_pos', 'gene','chr', 'end_pos')],
                            input_file[ ,c('gene', 'chr','pos', 'OrigDepth', 'Depth_30', "FwdSig", 'RevSig')],
                            by.y = c("gene","pos", "chr"), 
                            by.x = c('gene','trans_pos', 'chr'),
                            all.x = TRUE)
  get_text_exon_bed$Depth_30 = as.numeric(get_text_exon_bed$Depth_30)
  get_text_exon_bed$FwdSig = as.numeric(get_text_exon_bed$FwdSig)
  get_text_exon_bed$RevSig = as.numeric(get_text_exon_bed$RevSig)
  get_text_exon_bed$trans_pos = as.numeric(get_text_exon_bed$trans_pos)
  get_text_exon_bed = as.data.frame(get_text_exon_bed,stringsAsFactors=F)
  get_text_exon_bed$transcript_position = seq(1:length(get_text_exon_bed$gene))
  return(get_text_exon_bed)
}


##Makes graphs 
#Input file : text exon bed
# out put file : PDF written to the working directory 

graphs  = function(x){
  
  ablines = x[!duplicated(x["end_pos"]),]
  genename = unique(x$gene)
  plot(x$transcript_position, x$Depth_30, type="l", yaxt="n", col="black", ylab = "", xlab ="Transcript Position",main=genename)
  abline(h = c(50,100,200), lty = 'dashed', col = 'red', untf = T)
  axis(2, ylim=c(0,1))
  mtext(side =2, "Depth >Q30") 
  par(new=T)
  plot(x$transcript_position, x$FwdSig, type="l", ylim = c(0,1), yaxt="n", col= rgb(0.0,0.5,0.9, 0.3), ylab = "", xlab ="Transcript Position", lwd = 1)
  par(new= T)
  plot(x$transcript_position, x$RevSig, type = 'l', ylim = c(0,1), yaxt = 'n',col = rgb(1,0.8,0, 0.3), ylab = "", xlab ="Transcript Position", lwd =1)
  mtext(side = 4, "Variant signal")
  axis(side = 4, ylim = c(-1,1))
  legend("topright",  c("Depth", "Fwd_sig", "Rev_sig"), lwd= 1, col=c("black", rgb(0.0,0.5,0.9, 0.3), rgb(1,0.8,0, 0.3)), 
         lty = c(1,1,1), box.lty=0)
  abline(v = unlist(ablines$transcript_position), col = "darkgrey", lwd = 0.4, lty = 'dashed', untf = T)
  
}


##function to calculate regions in cosmic files
#input file : 
# output file: 
regions = function(cosmic_template){
  cosmic_file2 = cSplit(cosmic_template, "Mutation.Genome.Position", ":")
  cosmic_file2 = cSplit(cosmic_file2, "Mutation.Genome.Position_2", "-")
  return(cosmic_file2)
}


##function to calculate depth, regions and cosmic percentage 
#input file :
#output file : 
depth_cosmic_perctg = function(cosmic_region_template, text_exon_bed_template){
  colnames(cosmic_region_template) = c("Gene.name", "SNP", "chr", 
                                       "trans_pos", "Mutation.Genome.Position_2_2")
  try1 = text_exon_bed_template[text_exon_bed_template$trans_pos %in% 
                                  intersect(cosmic_region_template$Mutation.Genome.Position_2_2, 
                                            text_exon_bed_template$trans_pos),]
  
  cosmic_rep = as.data.frame(table(cosmic_region_template$Mutation.Genome.Position_2_2))
  try_2 = merge(try1, 
                cosmic_rep[, c('Var1', 'Freq')],
                by.x = c('trans_pos'),
                by.y = c('Var1'))
  n = length(cosmic_region_template$Gene.name)
  
  try_2$cosmic_prct = (try_2$Freq/n)*100
  
  colnames(try_2) = c("Start", "Exon", "Gene", "Chr", "End", "OrigDepth", "Depth30", "FwdSig", "RevSig", "TransPos", "Freq", "COSMIC_Pcrtg")
  try_2 = try_2[, c ("Gene", "Chr","Start", "End", "Exon", "Depth30", "Freq", "COSMIC_Pcrtg")]   
  return(try_2)
}

##added today 

##regions for bases with depth issues and cosmic cosmic percentage 

#  regions_depth_file = function(try_2_template){
#  regions_cos_depth = data.frame(seqToIntervals(try_2_template$trans_pos))
#  mrgerd_regions = merge(try_2_template, regions_cos_depth)
#  mrgerd_regions_filter = as.data.frame(filter(mrgerd_regions, trans_pos <= to & trans_pos >= from))
#  mrgerd_regions_filter_2 <- ddply(mrgerd_regions_filter
#                                   ,c('from','to', 'gene', 'chr', 'V5')
#                                   , summarize
#                                   , sum_cosmic_total = sum(cosmic_prct))
#  mergerd_regions <- within(mrgerd_regions_filter_2,  Region <- paste(from, to,  sep="-"))
#  merged_regions_2 = select(mergerd_regions, c("gene", "chr", "V5", "sum_cosmic_total", "Region"))
#  colnames(merged_regions_2)[1:5] = c("Gene", "Chr","Exon number","COSMIC percentage", "Region" )
#  merged_regions_2 = merged_regions_2[c("Region", "Gene", "Chr", "Exon number", "COSMIC percentage")]
#  return(merged_regions_2)  

#}





###########################***************************#############################

####Step1 :: Load tables 

generate_output_files = function(totgenelist,outputpdfname,ppdf_filepath,ppdf_input_large_file,dblocation_and_name){
  print(length(totgenelist$V1))
  pdf_name = outputpdfname
  
  pdf(pdf_name)
  
  
  pdf(file=paste0(outputpdfname,'.pdf'), width = 11, height = 8)
  par(mfrow= c(1,1))
  genedepthfile = data.frame()
  genedirectionalnoise = data.frame()
  perct_depth_cosmic_all = data.frame()
  perct_direct_noise_all = data.frame()
  gene_directional_noise_1 = data.frame() 
  percent_direct_noise_total = data.frame()
  depth_issues_cosmic_total = data.frame()
  depth_issues_total = data.frame()
  depth_cosmic_prctg_regions_all= data.frame()
  
  
  # Step :: working on genelist 
  for (i in 1:length(totgenelist$V1)) {
    
    
    print(1)  
    
    ##Step2 : query tables from SQL  
    gnnm=totgenelist$V1[i]
    text_file = small_data_1(gnnm,dblocation_and_name)
    exon_bed_file = small_data_2(gnnm,dblocation_and_name)
    
    #adding exon number 
    print(gnnm)
    print(nrow(exon_bed_file))
    exon_bed_file$V5 = (1:length(exon_bed_file$V1))
    
    cosmic_file3 = small_data_3(gnnm)
    #opening the cosmic file to separate chr, start, end into 3 columns 
    cosmic_file_2 = regions(cosmic_file3) 
    
    print(2)
    
    ##Step 3: work on exon file ( combine)
    
    ebf = exon_bed_file 
    df = data.frame()
    
    print(3-1)
    
    #Step 3-a :: fnding all the transcript postions using start and end using the bed file   
    df_2 = opening_exon_file(ebf)
    
    ##notes :df2 has all the postions mentioned in the exon bed file ( trans_pos, gene, chr, end_pos, exon)
    
    print(3-2)
    #step 3-b :: Merging (extracting data from the text file for the bed file)
    
    text_exon_bed = make_text_exon_bed(df_2, text_file)
    
    print(3-3)
    #step 3_c_1 :: creates df3 and
    
    df_3 = merge(unique(ebf[ ,c("V3", "V2","V1","V4","V5")])
                 , text_exon_bed[, c('trans_pos', "gene", "chr", 'transcript_position')],
                 by.x = c('V4', 'V1'),
                 by.y = c('gene', 'chr'),
                 all.y= TRUE)
    df_3 = as.data.frame(filter(df_3, trans_pos >= V2 & trans_pos <= V3))
    
    print(3-4)
    
    #Step3_c_2 :: adding exon number 
    
    text_exon_bed = merge(df_3[,c("V5","trans_pos")],text_exon_bed,by= c("trans_pos"))
    
    print(3-5)
    ##step 4 :: graphs
    
    graphs(text_exon_bed)
    
    print(4)
    
    ##Step 4 :: generating depth file 
    
    print(5)
    
    gene_depth_files_region = calculate_depth(text_exon_bed)
    if (class(gene_depth_files_region) == 'data.frame'){
      genedepthfile = rbind(genedepthfile, gene_depth_files_region ) 
    } else {
      genedepthfile  = " "
    }
    
    
    print(6)
    
    #Step 5_a :: Cosmic and region calculation 
    
    mergerd_regions_filter_final = depth_cosmic_perctg(cosmic_file_2, text_exon_bed)
    
    print("6.1.1")
    #depth_cosmic_prctg_regions = regions_depth_file(mergerd_regions_filter_final)
    depth_cosmic_prctg_regions_all= rbind(depth_cosmic_prctg_regions_all, mergerd_regions_filter_final)
    ##total bases with depth issues:
    print("6.1")
    if (class(gene_depth_files_region) == "data.frame"){
      perct_depth_1 = (nrow(gene_depth_files_region)/ nrow(text_exon_bed))*100
    } else {
      perct_depth_1  = 0
    }
    
    print("6.2")
    depth_issues_row = as.data.frame(cbind(Gene= text_exon_bed[2,3],metric='Percentage of bases with depth issues',
                                           value=perct_depth_1 ),stringsAsFactors=F)  
    depth_issues_total = rbind(depth_issues_total, depth_issues_row )
    print(10)
    
    
    
    
    
    ##total bases with depth issue and cosmic coverage 
    
    #if (class(mergerd_regions_filter_final) == "data.frame"){
    #  perct_depth_cosmic = (nrow(mergerd_regions_filter_final)/ nrow(text_exon_bed))*100
    #} else {
    #  perct_depth_cosmic  = 0
    # }
    
    # depth_issues_cosmic_row = as.data.frame(cbind(Gene= text_exon_bed[2,3],metric='No. of genes with depth issues and COSMIC coverage',
    #                                       value=perct_depth_cosmic),stringsAsFactors=F)  
    #  depth_issues_cosmic_total = rbind(depth_issues_cosmic_total, depth_issues_cosmic_row )
    print(10)
    
    
    
    print(7)
    #step5_b:: genetaing directional file 
    
    gene_directional_noise_merge = calculate_directional_noise(text_exon_bed)
    #gene_directional_noise_merge_2 = adding_regions(gene_directional_noise_merge)
    gene_directional_noise_1 = rbind(gene_directional_noise_1, gene_directional_noise_merge)
    
    
    print(8)
    #total bases with directional noise 
    
    
    perct_direct_noise = (nrow(gene_directional_noise_merge)/ nrow(text_exon_bed))*100
    direct_prctg_row = as.data.frame(cbind(Gene= text_exon_bed[2,3],metric='Percentage of bases with directional noise ',
                                           value=perct_direct_noise),stringsAsFactors=F)  
    
    percent_direct_noise_total = rbind(percent_direct_noise_total, direct_prctg_row )
    
    print(9)
    
    
    
    print(10)
    
    
    #return(depth_cosmic_prctg_regions_all)
    
  }
  
  ###LOOP ENDS HERE
  dev.off()
  o4 <- as.data.frame(cbind(Gene=" ",metric= " ",value= " "),stringsAsFactors=F)
  summary_table1 = rbind(percent_direct_noise_total,o4, o4, o4, depth_issues_total)
  
  
  
  print('a')
  ##total number of genes examined 
  total_genes_examined = length(totgenelist$V1)
  
  out1 = as.data.frame(cbind(Gene='All',metric='Total # Genes examined',
                             value=total_genes_examined),stringsAsFactors=F)  
  print('b')
  ## percentage of genes with depth issues
  
  genedepthfile_2 = genedepthfile
  
  genedepthfile_2 = genedepthfile_2[!apply(genedepthfile_2 == " .", 1, all),]
  
  freq = as.data.frame(table(genedepthfile_2$gene))
  print("done")
  total_gene_depth_issues = length(freq$Freq)
  print("done")
  
  percentage_depth_issues = (total_gene_depth_issues/length(totgenelist$V1))*100
  
  print("c")
  
  out2 = as.data.frame(cbind(Gene='All',metric='Percentage of genes with depth issues',
                             value = percentage_depth_issues),stringsAsFactors=F) 
  
  ##Total number of genes with directional noise 
  
  #directional noise 1 or greater 
  
  genes_direct_issues = as.data.frame(table(gene_directional_noise_1$gene))
  total_genes_with_direct_noise_1 = nrow(genes_direct_issues)/nrow(totgenelist)*100
  
  out3 = as.data.frame(cbind(Gene='All',metric='total_genes_with_direct_noise_1',
                             value = total_genes_with_direct_noise_1),stringsAsFactors=F) 
  
  #direction noise 2 or greater   
  print("d")
  
  total_genes_with_direct_noise_2 = total_bases_with_direct_noise(gene_directional_noise_1, 2)
  
  
  print("d1")
  out4 = as.data.frame(cbind(Gene='All',metric='total genes uni directional noise greater than 2',
                             value = total_genes_with_direct_noise_2),stringsAsFactors=F) 
  print("e")
  #direction noise greater than 3 
  total_genes_with_direct_noise_3 = total_bases_with_direct_noise(gene_directional_noise_1,3)
  out5 = as.data.frame(cbind(Gene='All',metric='total genes uni directional noise greater than 3',
                             value = total_genes_with_direct_noise_3),stringsAsFactors=F) 
  
  print("f")
  #direction noise greater than 4 and above 
  total_genes_with_direct_noise_4 = total_bases_with_direct_noise(gene_directional_noise_1,4)
  
  out6 = as.data.frame(cbind(Gene='All',metric='total genes uni-directional noise greater than 4',
                             value = total_genes_with_direct_noise_4),stringsAsFactors=F) 
  
  ##total genes with no issues 
  
  abs = rbind(genes_direct_issues,freq)
  abs_t = nrow(as.data.frame(table(abs$Var1)))
  no_issue_gene_prct = ((nrow(totgenelist) - abs_t)/ nrow(totgenelist))*100
  
  out7 = as.data.frame(cbind(Gene='All',metric='Percentage of genes with No issues',
                             value = no_issue_gene_prct),stringsAsFactors=F) 
  
  ##saving results table 
  print("g")
  outall= rbind(out1, out7, out2,out3, out4, out5, out6)
  print("h")
  
  
  
  
  write.table(summary_table1
              , paste0(ppdf_filepath,'/',"summary_table_1.txt"), sep="\t",
              quote =F, row.names=FALSE,
              col.names = colnames(summary_table1))
  
  
  write.table(outall
              , paste0(ppdf_filepath,'/', "summary_table_2.txt"), sep="\t",
              quote =F, row.names=FALSE,
              col.names = colnames(outall))
  
  write.table(depth_cosmic_prctg_regions_all
              , paste0(ppdf_filepath,'/',ppdf_input_large_file,".cosmic_depths.txt"), sep="\t",
              quote =F, row.names=FALSE,
              col.names = colnames(depth_cosmic_prctg_regions_all))
  
  write.table(gene_directional_noise_1, 
              paste0(ppdf_filepath,'/',ppdf_input_large_file,".direct_noise.txt"),
              sep="\t", quote = F,
              row.names=FALSE,
              col.names = colnames(gene_directional_noise_1))
  
  
}


################ FINAL FUNCTION
load_and_run_analysis = function(dblocation_and_name
                                 , filepath = getwd()
                                 , input_gene_file
                                 , output_pdf
                                 , input_large_file
                                 , input_bed_file
                                 , input_cosmic_file
                                 , reload_files = F
                                 , ...){
  
  print(paste('Working in',filepath,' with file ',input_large_file))
  
  if (reload_files){
    additional_largefiles = unlist(list(...))
    print(additional_largefiles)
    
    load_largefile_to_sql(input_large_file,dblocation_and_name,new_load = T)
    
    additional_largefiles = additional_largefiles[additional_largefiles!='NA']
    if (length(additional_largefiles) > 0){
      for (j in 1:length(additional_largefiles)){
        load_largefile_to_sql(additional_largefiles[j],dblocation_and_name,new_load = F)
      }
      
    }
    
    load_file_to_sql(input_bed_file
                     ,input_cosmic_file
                     ,dblocation_and_name)
    
  }
  gene_list = load_gene_list(paste0(filepath,'/',input_gene_file))
  
  textfiles =generate_output_files(gene_list,paste0(filepath,'/',output_pdf)
                                   ,filepath,input_large_file, dblocation_and_name)
  
  
}

if (length(input_args) == 7){
  load_and_run_analysis(dblocation_and_name = input_args[1]
                        , filepath = input_args[2]
                        , input_gene_file = input_args[3]
                        , output_pdf = input_args[4]
                        , input_large_file = input_files_arg5
                        , input_bed_file = args6
                        , input_cosmic_file = args7
                        , reload_files = T
                        , ... = input_files_arg5_ellipsis)
  
} else if (length(input_args) == 4)  {
  load_and_run_analysis(dblocation_and_name = input_args[1]
                        , filepath = input_args[2]
                        , input_gene_file = input_args[3]
                        , output_pdf = input_args[4]
                        , input_large_file = 'NVM'
                        , input_bed_file = 'NVM'
                        , input_cosmic_file = 'NVM'
                        , reload_files = F)
  
} else {
  print('Incorrect number of arguments  entered - please enter 4 args or 7 args')
}