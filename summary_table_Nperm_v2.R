######################################################################################
# CRG 
# Hana SUSAK
# date: 11/05/2016
#------------------------------------------------------------------------------------
# Make summary tables for N permutatios tests
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Calling the script: Rscript summary_table_Nperm_v2.R --h
#------------------------------------------------------------------------------------
######################################################################################

rm(list=ls())
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

#setwd('/no_backup/GD/projects/CLL/germ_analysis/')

#########################
parser$add_argument("-r", "--result_folder", type="character", help="Folder with resuls", metavar="folder path", nargs=1, required=TRUE)
parser$add_argument("-m", "--min_mut_cases", type="integer", help="Post Filter: number of minimum mutations in cases", metavar="number")
parser$add_argument("-p", "--positiv_sel", type="character", help="Post Filter: only when there are more affected cases then controls", metavar="T or F")

args <- parser$parse_args()

folder <- args$result_folder
# folder <- 'results/temp_hana/AF0.005/'
nn <- args$min_mut_cases
# nn <- 5
pos.sel <- as.logical(args$positiv_sel) 
if(!is.null(pos.sel)){
  if(is.na(pos.sel)){
    stop("-p (--positiv_sel) input parameter can be only T or F")
  } else {
    print("Have in mind you had to use same number of contorls and cases for this filter to have seanse")
  }
}

# pos.sel <- T
files <- list.files(folder, recursive=F, full.names=T)
files <- files[grepl('perm',files)]

ff<- list()
for (i in files){
  name <- unlist(strsplit(i,'/'))
  name <- name[length(name)]
  name <- unlist(strsplit(name,'_'))
  name <- name[1]
  print(name)
  ranking.i <- read.table(i, stringsAsFactors=F, sep='\t', header=T)
  
  # p vals stats
  p.val.cols <- which(grepl('p.val.overall_perm',colnames(ranking.i)))
  quantils.i <- apply(ranking.i[,p.val.cols],1, function(x)quantile(x, c(0.5,0.95), na.rm=T ))
  sig.p.n <- apply(ranking.i[,p.val.cols],1, function(x) sum(x <=0.05,na.rm=T))
  #not.na <- apply(ranking.i[,p.val.cols],1, function(x) sum(!is.na(x)))
  #df.p <- cbind.data.frame(as.data.frame(ranking.i$genes), as.data.frame(not.na), as.data.frame(sig.p.n), t(quantils.i))
  df.p <- cbind.data.frame(as.data.frame(ranking.i$genes), as.data.frame(sig.p.n), t(quantils.i))
  colnames(df.p) <- c("gene", "#sig_p", "quant_50_pval", "quant_95_pval")
  
  # adjusted p vals stats
  p.adj.val.cols <- which(grepl('p.adjust_perm',colnames(ranking.i)))
  quantils.i <- apply(ranking.i[,p.adj.val.cols],1, function(x) quantile(x, c(0.5,0.95), na.rm=T ))
  sig.p.adj.n <- apply(ranking.i[,p.adj.val.cols],1, function(x) sum(x <=0.1,na.rm=T))
  df.adj <- cbind.data.frame(as.data.frame(ranking.i$genes),  as.data.frame(sig.p.adj.n), t(quantils.i))
  colnames(df.adj) <- c("gene", "#sig_p_adj","quant_50_pval_adj",  "quant_95_pval_adj")   
  
  # controls coutns
  contrl.cols <- which(grepl('num.controls_perm',colnames(ranking.i)))
  quantils.i <- apply(ranking.i[,contrl.cols],1, function(x)quantile(x, c(0.5,1), na.rm=T ))
  df.controls <- cbind.data.frame(as.data.frame(ranking.i$genes),  t(quantils.i))
  #colnames(df.controls) <- c("gene", "quant_0_cont", "quant_25_cont", "quant_50_cont", "quant_75_cont", "quant_100_cont")
  colnames(df.controls) <- c("gene", "quant_50_cont",  "quant_100_cont")
  
  # cases coutns
  cases.cols <- which(grepl('num.cases_perm',colnames(ranking.i)))
  quantils.i <- apply(ranking.i[,cases.cols],1, function(x) quantile(x, c(0.5,1), na.rm=T ))
  df.cases <- cbind.data.frame(as.data.frame(ranking.i$genes),  t(quantils.i))
  #colnames(df.cases) <- c("gene", "quant_0_cases", "quant_25_cases", "quant_50_cases", "quant_75_cases", "quant_100_cases")
  colnames(df.cases) <- c("gene", "quant_50_cases",  "quant_100_cases")
  
  df.fin <- cbind.data.frame(df.p, df.adj[,-1], df.cases[,-1], df.controls[-1])
  df.fin <- df.fin[order( df.fin$`#sig_p_adj`, df.fin$`#sig_p`, -df.fin$quant_95_pval_adj,-df.fin$quant_95_pval, df.fin$quant_50_cases, -df.fin$quant_50_cont ,decreasing=T),]
  
  if(!is.null(pos.sel)){
    df.fin <- df.fin[df.fin$quant_50_cases > df.fin$quant_50_cont,]
  }
  if(!is.null(nn)){
    df.fin <- df.fin[df.fin$quant_50_cases > nn,]
  }
  
  ff[[name]] <- df.fin
}


rm(df.fin)
rm(df.adj)
rm(df.cases)
rm(df.controls)
rm(df.p)
rm(quantils.i)
rm(ranking.i)
rm(name)
rm(i)
rm(cases.cols)
rm(contrl.cols)
#rm(not.na)
rm(p.val.cols)
rm(p.adj.val.cols)
rm(sig.p.n)
rm(sig.p.adj.n)


for (f in 1:length(ff)){
  name <- unlist(strsplit(files[f], '/'))
  name <- name[length(name)]
  name2 <- names(ff)[f]
  print(name2)
  
  if(!is.null(pos.sel)){
    write.table(ff[[f]], paste0(name,'_post_Filtered_summary.txt'), sep='\t', quote=F, col.names=T, row.names=F )
  } else {
    write.table(ff[[f]], paste0(name,'_summary.txt'), sep='\t', quote=F, col.names=T, row.names=F )
  }
}


rm(name)
