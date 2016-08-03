######################################################################################
# CRG 
# Hana SUSAK
# date: 01/02/2016
#------------------------------------------------------------------------------------
# SKAT-O, MiST and KBAC test, with or without permutations for controls
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Calling the script: Rscript germline_risk_test_v2.R --h
#------------------------------------------------------------------------------------
######################################################################################

rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MiST))
suppressPackageStartupMessages(library(SKAT))
suppressPackageStartupMessages(library(KBAC,lib.loc="/software/xe/el6.3/R_libs"))
suppressPackageStartupMessages(library(parallel))


parser <- ArgumentParser()

parser$add_argument("-m", "--mutation.file", type="character", help="input snp/indel file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-d", "--sample_desc_file", type="character", help="samples pca and description file (first column sample ID)", metavar="file", nargs=1, required=TRUE)
#parser$add_argument("-p", "--pca_col_num", type="integer", help="Number of PCA components in samples info file (starting from second column)", metavar="number", required=TRUE)
#parser$add_argument("-p", "--pca_results", type="character", help="file with PCA results", metavar="file", nargs=1,required=TRUE)
parser$add_argument("-f", "--allele_freq", type="character", help="Allele frequency threshold", metavar="file", nargs=1,  default="0.01")
parser$add_argument("-t", "--genes_list_file", type="character", help="Genes to test file", metavar="file", nargs=1)
#parser$add_argument("-c", "--correct_sample_names", type="character", help="file with corrected names", metavar="file", nargs=1)
parser$add_argument("-o", "--out_folder", type="character", help="Folder to set as working directory", metavar="directory", nargs=1)


read.from <- as.character(Sys.getenv("READ_FROM", 'stdin')) 
# read.from <- ''
cat('Standard input is:',read.from,'\n')
if (read.from=='bash'){
   mc <- as.integer(Sys.getenv("R_MAX_MC_CORES", 1L))    
   seed <- as.integer(Sys.getenv("R_SEED", 160185L)) 
   set.seed(seed)
   MIST_log <- as.logical(Sys.getenv("R_MIST", 'T')) 
   KBAC_log <- as.logical(Sys.getenv("R_KBAC", 'T')) 
   SKATO_log <- as.logical(Sys.getenv("R_SKATO", 'T')) 
} else if (read.from == 'stdin'){ 
   cat("\nEnter a number of cores available: ")
   mc <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nEnter a seed for random slecting of control samples: ")
   seed <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   set.seed(seed)
   cat("\nDo you want to run MiST (T for yes, F for no): ")
   MIST_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nDo you want to run KBAC (T for yes, F for no): ")
   KBAC_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nDo you want to run SKAT-O (T for yes, F for no): ")
   SKATO_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
   
} else {
   mc <- as.integer(1L) 
   seed <- as.integer(160185L) 
   set.seed(seed)
   MIST_log <- T
   KBAC_log <- T
   SKATO_log <- T
}



cat("\n Number of cores useing", mc, " .....\n" )
cat("\n Seed set to", seed, "..... \n " )
cat(paste0("Tests to be performed: ", paste(c('MiST','KBAC','SKAT-O')[c(MIST_log,KBAC_log,SKATO_log)], collapse=", ")), "\n")

######################################################################################
# read arguments
######################################################################################

## get command line arguments
args <- parser$parse_args()
mutation.file <- args$mutation.file
samples.info.file <- args$sample_desc_file
af.thr <- as.numeric(args$allele_freq)
genes.file <- args$genes_list_file
out.folder <- args$out_folder
#pca.cols <- args$pca_col_num
#pca.file <- args$pca_results
#sample.correct.file.name <- args$correct_sample_names 

if(!is.null(out.folder)){
   cat('Setting to out folder :',out.folder, '\n')
   setwd(out.folder)
}

# mutation.file <- '/no_backup/GD/projects/CLL/germ_analysis/data/QC_output/mutations_filtered_samples_fix.txt'
# samples.info.file <- '/no_backup/GD/projects/CLL/germ_analysis/data/QC_output/samples_info_pca.txt'
# pca.cols <- 20
# af.thr <- 0.005
# genes.file <- '/no_backup/GD/projects/CLL/germ_analysis/results/AF0.005/genes_top_43.txt'
cat('Mutations file is:',mutation.file,'\n')
cat('Samples info file is:',samples.info.file,'\n')

#creat log file
logFile <-  paste0('risk_analysis_log_file_',format(Sys.time(), format = "%Y-%j-%H%M%S") , '.txt')
cat(timestamp(), file=logFile, append=FALSE, sep = "\n")

cat(paste0("Number of cores: ", mc), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Seed: ", seed), file=logFile, append=TRUE, sep = "\n")

cat(paste0("Tests to be performed: ", paste(c('MiST','KBAC','SKAT-O')[c(MIST_log,KBAC_log,SKATO_log)], collapse=", ")),
    file=logFile, append=TRUE, sep = "\n")

cat(paste0("Mutations file: ", mutation.file), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Samples info file: ", samples.info.file), file=logFile, append=TRUE, sep = "\n")


######################################################################################
# functions 
######################################################################################
if(MIST_log) {
   doit.mist <- function(y.bin, X, G, Z.func, maf){
      #print(paste('Gene ',gene, 'is tasted'))
      test.ex <- logit.weight.test(y.bin,X,G,Z.func,maf)
      return(test.ex)
   } 
   
   test.mist.logit <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, pi){   
      write(paste('starting to test gene: ',gene), file=paste('mist_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( '#Gene' == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#SpanishFreq')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         
         mut.total <- ncol(G)
         mut.cases <- sum(colSums(as.matrix(G[cases.samp ,])) !=0)
         mut.controls <- sum(colSums(as.matrix(G[controls.samp,])) !=0)
         num.cases <- sum(rowSums(as.matrix(G[cases.samp ,])) !=0)
         num.controls <- sum(rowSums(as.matrix(G[controls.samp,])) !=0)
         
         maf.df <- pop.mut[ colnames(G),]         
         maf <- maf.df$'#SpanishFreq' + 0.0001
         names(maf)  <- maf.df$'#snpID'
         rm(maf.df)
         
         y.bin <- as.numeric(rownames(G) %in% cases.samp)
         
         Z.func <- data.frame(intercept=1, cadd=mutations[ colnames(G),'#Cadd2'])
         
         Z.func[is.na(Z.func$cadd),'cadd'] <- 15
         rownames(Z.func) <- colnames(G)
         
         test.ex <- (try( doit.mist(y.bin, Xi, G, Z.func, maf), TRUE))
         
         if(length(test.ex) == 5){
            p.val.add <-  as.numeric(test.ex$p.value.overall)
            p.val.tau <-  as.numeric(test.ex$p.value.S.tau)
            p.val.pi <- as.numeric(test.ex$p.value.S.pi)
         } else {
            p.val.add <-  NA
            p.val.tau <- NA
            p.val.pi <- NA
         } 
      } else {
         mut.total <- 0
         mut.cases <- 0
         mut.controls <- 0
         num.cases <-0
         num.controls <- 0
         p.val.add <-  NA
         p.val.tau <- NA
         p.val.pi <- NA
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls, p.val.pi=p.val.pi,  p.val.tau=p.val.tau,
                 p.val.overall=p.val.add))     ) 
   }
   
   perm.mc.mist <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per, pi) {
      write('starting for this permutation!', file=paste('mist_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.mist.logit(x, mutations.per[mutations.per$'#Gene'==x,],
                                                          cases.samples,controls.samples.per,
                                                          pop.mut.df.per, X.per, pi))
      
      write('Done with this permutation!', file=paste('mist_perm',pi, sep='_'), sep = "\t", append=T)
      
      df.temp <- as.data.frame(t(df))
      
      df.temp <- data.frame(df.temp, p.adjust = p.adjust(df.temp$p.val.overall, "BH", n=length(genes.all)))
      
      colnames(df.temp) <- paste(colnames(df.temp), '_perm',pi, sep='')
      
      write('going out of perm.mc.mist function!', file=paste('mist_perm',pi, sep='_'), sep = "\t", append=T)
      
      return(df.temp)
      
   }
   
}

####

if(SKATO_log) {
   doit.skat <- function(y.bin, X, G, maf){
      
      obj.s <- SKAT_Null_Model(y.bin ~ X, out_type="D",Adjustment=FALSE)
      weights <- Get_Logistic_Weights_MAF(maf)   
      test.ex <- SKATBinary(G,obj.s,kernel="linear.weighted",weights=weights,method="optimal.adj")
      
      return(test.ex)
   } 
   
   test.skat.binary <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, pi){   
      write(paste('starting to test gene: ',gene), file=paste('skat0_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( '#Gene' == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#SpanishFreq')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         
         mut.total <- ncol(G)
         mut.cases <- sum(colSums(as.matrix(G[cases.samp ,])) !=0)
         mut.controls <- sum(colSums(as.matrix(G[controls.samp,])) !=0)
         num.cases <- sum(rowSums(as.matrix(G[cases.samp ,])) !=0)
         num.controls <- sum(rowSums(as.matrix(G[controls.samp,])) !=0)
         
         maf.df <- pop.mut[ colnames(G),]         
         maf <- maf.df$'#SpanishFreq' + 0.001
         names(maf)  <- maf.df$'#snpID'
         rm(maf.df)
         
         y.bin <- as.numeric(rownames(G) %in% cases.samp)
         
         
         #if( (mut.cases>0 & mut.controls>0) | (mut.cases==0 & num.controls>2) | (mut.controls==0 & num.cases>2)){
         test.ex <- (try( doit.skat(y.bin, as.matrix(Xi), G, maf), TRUE))
         #} else {
         #   test.ex <- NA
         #}
         
         if(length(test.ex) >= 1){
            p.val <-  test.ex$p.value  
         } else {
            p.val <-  NA
         }
      } else {
         mut.total <- 0
         mut.cases <- 0
         mut.controls <- 0
         num.cases <-0
         num.controls <- 0
         p.val <-  NA
         
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls, 
                 p.val.overall=p.val))     ) 
   }
   
   perm.mc.skato <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per, pi) {
      write('starting for this permutation!', file=paste('skat0_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.skat.binary(x, mutations.per[mutations.per$'#Gene'==x,],
                                                           cases.samples,controls.samples.per,
                                                           pop.mut.df.per, X.per, pi))
      
      write('Done with this permutation!', file=paste('skat0_perm',pi, sep='_'), sep = "\t", append=T)
      
      df.temp <- as.data.frame(t(df))
      
      df.temp <- data.frame(df.temp, p.adjust = p.adjust(df.temp$p.val.overall, "BH", n=length(genes.all)))
      
      colnames(df.temp) <- paste(colnames(df.temp), '_perm',pi, sep='')
      
      write('going out of perm.mc.skato function!', file=paste('skat0_perm',pi, sep='_'), sep = "\t", append=T)
      
      return(df.temp)
      
   }
   
}
   
####

if(KBAC_log) {
   doit.kbac <- function(y.bin,G ){
      
      casectrl.dat <- as.matrix(cbind(y.bin,G))
      
      alpha <- 0.05
      num.permutation <- 3000
      quiet <- T
      alternative <- 1
      maf.upper.bound <- 1
      
      test.ex <- KbacTest(casectrl.dat, alpha, num.permutation, quiet, maf.upper.bound, alternative)
      
      return(test.ex)
   } 
   
   test.kbac <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, pi){   
      write(paste('starting to test gene: ',gene), file=paste('kbac_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( '#Gene' == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#SpanishFreq')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         
         mut.total <- ncol(G)
         mut.cases <- sum(colSums(as.matrix(G[cases.samp ,])) !=0)
         mut.controls <- sum(colSums(as.matrix(G[controls.samp,])) !=0)
         num.cases <- sum(rowSums(as.matrix(G[cases.samp ,])) !=0)
         num.controls <- sum(rowSums(as.matrix(G[controls.samp,])) !=0)
         
         maf.df <- pop.mut[ colnames(G),]         
         maf <- maf.df$'#SpanishFreq' + 0.001
         names(maf)  <- maf.df$'#snpID'
         rm(maf.df)
         
         y.bin <- as.numeric(rownames(G) %in% cases.samp)
         
         
         test.ex <- (try( doit.kbac(y.bin, G), TRUE))
         
         
         if(length(test.ex) > 0){
            p.val <-  test.ex         
         } else {
            p.val <-  NA
         } 
      } else {
         mut.total <- 0
         mut.cases <- 0
         mut.controls <- 0
         num.cases <-0
         num.controls <- 0
         p.val <-  NA
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls,
                 p.val.overall=p.val))) 
   }
   
   perm.mc.kbac <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per, pi) {
      write('starting for this permutation!', file=paste('kbac_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.kbac(x, mutations.per[mutations.per$'#Gene'==x,],
                                                    cases.samples,controls.samples.per,
                                                    pop.mut.df.per, X.per, pi))
      
      write('Done with this permutation!', file=paste('perm',pi, sep='_'), sep = "\t", append=T)
      
      df.temp <- as.data.frame(t(df))
      
      df.temp <- data.frame(df.temp, p.adjust = p.adjust(df.temp$p.val.overall, "BH", n=length(genes.all)))
      
      colnames(df.temp) <- paste(colnames(df.temp), '_perm',pi, sep='')
      
      write('going out of perm.mc.kbac function!', file=paste('kbac_perm',pi, sep='_'), sep = "\t", append=T)
      
      return(df.temp)
      
   }
   
}


######################################################################################

######################################################################################
# read input
######################################################################################
# mutations file
mutations.df <- read.table(file=mutation.file, header=T, comment.char="", quote="", sep='\t', #nrows=10000,
                           stringsAsFactors=F, check.names=F)

# - get all samples in mutation.df
annot.cols <-which(grepl('#',colnames(mutations.df)))
samples.cols <- which(!grepl('#',colnames(mutations.df)))
sample.start <- which('#snpID'==colnames(mutations.df))+1

mutations.df <- mutations.df[rowSums(mutations.df[,samples.cols])!=0 ,]
cat("\nIn input file there are", nrow(mutations.df), "unique snps/indels  \n ....." )

# samples info file
samples.info <- read.table(file=samples.info.file, header=T, sep="\t", quote="")
cat("\nThere is info for", nrow(samples.info), "samples. \n First column should be unique sample ID, as in multi-sample file \n ....." )

colnames(samples.info)[1] <- 'samples.my'

all.samples <- unique(as.character(colnames(mutations.df[,samples.cols])))
if ( all(sapply(all.samples, function(x) x %in% samples.info$samples.my)) ){
   samples.info <- samples.info[samples.info$samples.my %in% all.samples,]  
} else {
   stop("There must be info for all samples from mutation file. Some samples are mismached or missing.")
}

########################
# X matrix infos

if (read.from=='bash'){
   pca.cols <- as.integer(Sys.getenv("R_COV_MAT_COLS", 2L))    
   pca.x <- samples.info[,1:(pca.cols+1)]
   rownames(pca.x) <- pca.x$samples.my
   
} else if (read.from == 'stdin'){ 
   cat("\nThere are columns: ", 
       paste('\n',1:(ncol(samples.info)-1),'.',colnames(samples.info)[-1], sep=''),
       "\nWhich one you want to keep included in X (cov mat) matrix: ")
   x <- scan(read.from,what=integer(),nmax=ncol(samples.info)-1,quiet=TRUE, nlines=1)
   x <- c(1,x+1)
   pca.x <- samples.info[,x]
   rownames(pca.x) <- pca.x$samples.my
   rm(x)
} else {
   pca.cols <- 0    
   pca.x <- samples.info[,2:(pca.cols+1)]   
}

pca.x <- pca.x[pca.x$samples.my %in% all.samples , ]

########################
# AF threshold 
print(paste('Allele frequency threshold is: ',af.thr))  
cat(paste0("Allele frequency threshold: ", af.thr), file=logFile, append=TRUE, sep = "\n")

########################
# read genes to test
if(!is.null(genes.file)){
   print(paste('Genes list file name: ',genes.file))  
   genes.all <- read.table(file=genes.file, header=F, comment.char="", quote="", sep='\t', stringsAsFactors=F, check.names=F)  
   genes.all <- genes.all[!genes.all %in% c('Gene','gene','Gene_Symbol', '#Gene')] # in case of a header
   if (is.data.frame(genes.all)) {
      genes.all <- unique(as.character(genes.all[,1]))
   }
} else {
   genes.all <- unique(as.character(mutations.df$'#Gene'))  
}   
genes.all <- as.character(sort(genes.all))

cat("\nWe will test for : ", length(genes.all), "genes\n ...")
cat("\nWe will use in total: ", nrow(pca.x), "patients (some cases, some controls, some as af estimators )\n ...")
cat(paste0("Number of Genes to be tested: ", length(genes.all)), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Number of Samples to be tested (case+control): ",  nrow(pca.x)), file=logFile, append=TRUE, sep = "\n")

######################################################################################
# Samples for calculating local population AF
######################################################################################

if (read.from=='bash'){
   proj.col <- as.character(Sys.getenv("R_PROJECT_COL", 'DB_Project'))    
   cases.name <- as.character(Sys.getenv("R_CASES_NAME", 'CLL'))
   
   x <- which(proj.col == colnames(samples.info))
   options <- as.character(unique(samples.info[,x]))
   cases.names.num <- which(cases.name == options)
   
} else if (read.from == 'stdin'){ 
   cat("\n\nWhich column represents cases/controls (Projects)):", 
       paste('\n',1:ncol(samples.info),'.',colnames(samples.info), sep=''))
   cat("\nEnter column number: ")
   x <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nColumn(s) you chose is:",paste(colnames(samples.info)[x], collapse=', ' ), "\n.....")
   options <- as.character(unique(samples.info[,x]))
   options[is.na(options)] <- 'NA_Missing_info'
   cat("\n\nIn column",colnames(samples.info)[x],"there are folowing fields:")
   counts.tb <- table(samples.info[,x], useNA='always')
   names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
   counts.tb <- counts.tb[options]
   cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' samples)','\n', sep=''))
   cat("\nEnter  numbers which are considered as cases: ")
   cases.names.num <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
   cases.names.num <- unique(cases.names.num)
   cases.names.num <- cases.names.num[cases.names.num <= length(options)]
   cat("\nSamples whos", colnames(samples.info)[x],  "is/are", paste(options[cases.names.num], collapse=', '),
       "are considered as cases. \n.....") 
   #rm(x)   
   
} else {
   proj.col <-'DB_Project'    
   cases.name <-  'CLL'
   
   x <- which(proj.col == colnames(samples.info))
   options <- as.character(unique(samples.info[,x]))
   cases.names.num <- which(cases.name == options)
}

cat(paste0("Cases defined as: ", paste(options[cases.names.num], collapse=', ')), file=logFile, append=TRUE, sep = "\n")

# - split cases, controls and AF-set
cases.samples <- as.character(samples.info[samples.info[,x] %in% options[cases.names.num], 'samples.my' ])
controls.samples.all <- as.character(samples.info[!samples.info[,x] %in% options[cases.names.num], 'samples.my' ])

# chooce number of permutations for selecting cases and AF samples
if (read.from=='bash'){
   permut <- as.integer(Sys.getenv("R_PERMUTATIONS", 100L))    
   if ( permut < 1){
      stop("You should have number of permutations larger then 0")  
   }
   cadd.trh <- as.integer(Sys.getenv("R_CADD_TRH", 10L))    
   
} else if (read.from == 'stdin'){ 
   cat("\nHow menay controls permutations you want: ")
   permut <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   if ( permut < 1){
      stop("You should enter integer number larger then 0")  
   }
   
   cat("\nWhat CADD2 threshold you want to use as filter (put 0 if you want to keep all the mutations):")
   cadd.trh <- scan(read.from,what=numeric(),nmax=1,quiet=TRUE, nlines=1)
   
} else {
   permut <- 100
   cadd.trh <- 0
}
cat("\nCADD score as threshold used : ", cadd.trh, "\n ...")
cat(paste0("CADD score  threshold: ", cadd.trh), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Number of permutations: ", permut), file=logFile, append=TRUE, sep = "\n")

df.all <- data.frame()

same.perm <- FALSE

af.samples <- list()
controls.samples <- list()

for (p in 1:ceiling(permut/mc)){
   pp <- ((p-1)*mc+1):(p*mc)
   cat("\nPermutaion number:",pp,"!\n.....")
   # split controls- for local AF estimation
   if (p==1 | !same.perm ){
      if (read.from=='bash'){
         controls.num <- as.integer(Sys.getenv("R_CONTROLS_NUM", length(cases.samples)))  
         cat("\n",controls.num, "samples are used as contorols ...\n")
         
         answ <- as.character(Sys.getenv("R_PERM_ANSW", 'y'))    
      } else if (read.from == 'stdin'){ 
         cat("\nYou have",length(cases.samples),"cases samples, which leave you",length(controls.samples.all),
             "controls. \nHow many samples you want to use for controls (rest is use for estimating local polulation AF):")
         controls.num <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
         cat("\nYou want to use this number of controls for all permutations (answer should be y or n):")
         answ <- scan(read.from,what=character(),n=1,quiet=TRUE)
      } else {
         controls.num <- length(cases.samples)
         answ <- 'y'
      }
      cat(paste0("Number of controls used: ", controls.num), file=logFile, append=TRUE, sep = "\n")
      
      if(answ=='y'){
         same.perm <-  TRUE
      }       
   }
   

   if (controls.num < length(controls.samples.all)){
     
      for (ii in pp){
         controls.samples[[ii]] <- sort(sample(controls.samples.all,controls.num))
         af.samples[[ii]] <- controls.samples.all[!controls.samples.all %in% controls.samples[[ii]]] 
      }
      
   } else if (controls.num == length(controls.samples.all)) {
      controls.samples <- list()
      af.samples <- list()
      for (ii in pp){
         controls.samples[[ii]] <- controls.samples.all
         af.samples[[ii]] <- controls.samples.all 
      }   
   } else {
      stop('You have to choose number of cases less then available number of total cases')
   }

   pop.mut.df <- list()   
   for(ii in pp){
      # - calculate local AF from mutations.df (set of all snps maybe before)
      pop.cols <- which(colnames(mutations.df) %in% af.samples[[ii]])
      pop.mut <- mutations.df[,c( annot.cols , pop.cols)]
      rownames(pop.mut) <- pop.mut$`#snpID`
      freq <- rowSums(pop.mut[,-annot.cols])/(2*length(af.samples[[ii]]))
      # pop.mut$`#Freq` <- apply(cbind(freq,(1-freq)),1,min) # corrected minor af
      pop.mut$`#Freq` <- sapply(freq, function(x) min(x, 1-x) ) # corrected minor af      
      
      pop.samp.annot <- which(grepl('#',colnames(pop.mut)))
      pop.samp.cols <- which(!grepl('#',colnames(pop.mut)))
      pop.mut <- pop.mut[,c(pop.samp.annot, pop.samp.cols)]
      
      pop.mut.df[[ii]] <- pop.mut[,c('#snpID', '#Freq')]
      colnames(pop.mut.df[[ii]]) <- c('#snpID', '#SpanishFreq')
      #samples.start.pop <- sample.start+1
      rm(pop.cols)
      rm(pop.mut) 
   }   
   
   ########################
   #  fileter mutations 
   ########################  
   
   mutations.dam <- list()
   for(ii in pp){
      # 1KG eur frequency below af.thr
      mutations.df.fil <- mutations.df[mutations.df$'#Eur1000GenomesFrequency' < af.thr,]
      
      # EVS eur frequency below af.thr
      mutations.df.fil <- mutations.df.fil[mutations.df.fil$'#EurEVSFrequency' < af.thr,]
      
      # Spanish frequency below af.thr
      mutations.df.fil <- merge(x=mutations.df.fil, y=pop.mut.df[[ii]], by='#snpID', all.x=T)
      
      mut.fil.samp.annot <- which(grepl('#',colnames(mutations.df.fil)))
      mut.fil.samp.cols <- which(!grepl('#',colnames(mutations.df.fil)))
      mutations.df.fil <- mutations.df.fil[,c(mut.fil.samp.annot,mut.fil.samp.cols)]      
      mut.fil.samp.annot <- which(grepl('#',colnames(mutations.df.fil)))
      mut.fil.samp.cols <- which(!grepl('#',colnames(mutations.df.fil)))
      mut.fil.samp.cols.cases <- which(colnames(mutations.df.fil) %in% cases.samples)
      mut.fil.samp.cols.contorls <- which(colnames(mutations.df.fil) %in% controls.samples[[ii]])
      #mutations.df.fil <- mutations.df.fil[,c(1:(sample.start-1), ncol(mutations.df.fil),sample.start:(ncol(mutations.df.fil)-1))]
      sample.start.fil <- which('#SpanishFreq'==colnames(mutations.df.fil)) +1 
      
      mutations.df.fil <- mutations.df.fil[mutations.df.fil$`#SpanishFreq` < af.thr,]
      
      # nonsilent mutations
      mutations.df.fil <- mutations.df.fil[mutations.df.fil$`#ExonicFunction` != 'synonymous SNV' | 
                                              is.na(mutations.df.fil$`#ExonicFunction`), ]
      
      # correct indels cadd score if missign
      #mutations.df.fil[is.na(mutations.df.fil$'#Cadd2') & mutations.df.fil$'#type' == 'indel', '#Cadd2'] <- 40
      
      # only high CADD score
      mutations.df.fil <- mutations.df.fil[mutations.df.fil$'#Cadd2' >= cadd.trh | is.na(mutations.df.fil$'#Cadd2') | mutations.df.fil$'#type' == 'indel', ]
      
      #exclude samples for estimating local af
      if (controls.num < length(controls.samples.all)){
         mutations.df.fil <- mutations.df.fil[,!colnames(mutations.df.fil) %in%  af.samples[[ii]] ]
      }
      
      # take only rows which have at least one mutationin in cases and controls
      mut.fil.samp.cols <- which(!grepl('#',colnames(mutations.df.fil)))
      ind <- rowSums(mutations.df.fil[,mut.fil.samp.cols], na.rm=T)
      rownames(mutations.df.fil) <- mutations.df.fil$'#snpID'
      mutations.dam[[ii]] <- mutations.df.fil[ind!=0,]
      cat('In permutation', ii, 'there is', nrow(mutations.dam[[ii]]) ,' snps/indels ...\n' )
      rm(mutations.df.fil)
      
   } 
   
   # X matrix, same for all genes
   X <- list()
   for(ii in pp){
      Xi <- pca.x[c(cases.samples,controls.samples[[ii]]),]
      Xi$intercept <- 1
      Xi<- Xi[,c(ncol(Xi),2:(ncol(Xi)-1))] 
      X[[ii]] <- Xi
      rm(Xi)
   }
   
   
   ##############################################################################################
   # MIST
   ##############################################################################################   
   if(MIST_log){
      cat("\n Starting test mist!\n.....")
      df.perm.mist <- mclapply(pp, function(pi) perm.mc.mist(genes.all, mutations.dam[[pi]],
                                                             cases.samples, controls.samples[[pi]],
                                                             pop.mut.df[[pi]], X[[pi]], pi), mc.cores=mc)
      
      
      
      
      if (p == 1 ){
         df.all.mist <- Reduce(merge, lapply(df.perm.mist, function(x) data.frame(x, genes = row.names(x))))      
         
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.mist) ))
         if (length(cc)>1 ){
            df.all.mist2 <- df.all.mist[order(rowMeans(df.all.mist[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.mist2 <- df.all.mist[order(df.all.mist[,cc], decreasing=F),]
         }
         
         #head(df.all.mist2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         #out.file <- "MiST"
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('MiST_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.mist2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.mist2)
         
      } else {
         df.temp <- Reduce(merge, lapply(df.perm.mist, function(x) data.frame(x, genes = row.names(x))))  
         df.all.mist <- merge(x=df.all.mist, y=df.temp, by="genes", all=T)   
         rm(df.temp)
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.mist) ))
         if (length(cc)>1 ){
            df.all.mist2 <- df.all.mist[order(rowMeans(df.all.mist[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.mist2 <- df.all.mist[order(df.all.mist[,cc], decreasing=F),]
         }
         
         #head(df.all.mist2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('MiST_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.mist2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.mist2)
      }
      
   }

   ##############################################################################################
   # SKAT-O
   ##############################################################################################
   if(SKATO_log){
      cat("\n Starting skato!\n.....")
      df.perm.skato <- mclapply(pp, function(pi) perm.mc.skato(genes.all, mutations.dam[[pi]],
                                                               cases.samples, controls.samples[[pi]],
                                                               pop.mut.df[[pi]], X[[pi]], pi), mc.cores=mc)
      
      if (p == 1 ){
         df.all.skato <- Reduce(merge, lapply(df.perm.skato, function(x) data.frame(x, genes = row.names(x))))      
         
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.skato) ))
         if (length(cc)>1 ){
            df.all.skato2 <- df.all.skato[order(rowMeans(df.all.skato[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.skato2 <- df.all.skato[order(df.all.skato[,cc], decreasing=F),]
         }
         
         #head(df.all.skato2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         #out.file <- "SKAT_O"
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('SKAT_O_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.skato2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.skato2)
         
      } else {
         df.temp <- Reduce(merge, lapply(df.perm.skato, function(x) data.frame(x, genes = row.names(x))))  
         df.all.skato <- merge(x=df.all.skato, y=df.temp, by="genes", all=T)   
         rm(df.temp)
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.skato) ))
         if (length(cc)>1 ){
            df.all.skato2 <- df.all.skato[order(rowMeans(df.all.skato[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.skato2 <- df.all.skato[order(df.all.skato[,cc], decreasing=F),]
         }
         
         #head(df.all.skato2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('SKAT_O_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.skato2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.skato2)
      }
   }
   
  
   
   ##############################################################################################
   # KBAC
   ##############################################################################################
   if(KBAC_log) {
      cat("\n Starting kbac!\n.....")
      df.perm.kbac <- mclapply(pp, function(pi) perm.mc.kbac(genes.all, mutations.dam[[pi]],
                                                             cases.samples, controls.samples[[pi]],
                                                             pop.mut.df[[pi]], X[[pi]], pi), mc.cores=mc)
      
      if (p == 1 ){
         df.all.kbac <- Reduce(merge, lapply(df.perm.kbac, function(x) data.frame(x, genes = row.names(x))))      
         
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.kbac) ))
         if (length(cc)>1 ){
            df.all.kbac2 <- df.all.kbac[order(rowMeans(df.all.kbac[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.kbac2 <- df.all.kbac[order(df.all.kbac[,cc], decreasing=F),]
         }
         
         #head(df.all.kbac2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         #out.file <- "SKAT_O"
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('KBAC_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.kbac2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.kbac2)
         
      } else {
         df.temp <- Reduce(merge, lapply(df.perm.kbac, function(x) data.frame(x, genes = row.names(x))))  
         df.all.kbac <- merge(x=df.all.kbac, y=df.temp, by="genes", all=T)   
         rm(df.temp)
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.kbac) ))
         if (length(cc)>1 ){
            df.all.kbac2 <- df.all.kbac[order(rowMeans(df.all.kbac[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.kbac2 <- df.all.kbac[order(df.all.kbac[,cc], decreasing=F),]
         }
         
         #head(df.all.kbac2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('KBAC_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.kbac2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.kbac2)
      }
   }
   
   
}

save(af.samples, file='af.samples.Rdata')
save(controls.samples, file='controls.samples.Rdata')

cat(paste0("In file af.samples.Rdata are saved (R object) af samples IDs used for AF estimation"), file=logFile, append=TRUE, sep = "\n")
cat(paste0("In file controls.samples.Rdata are saved (R object) af samples IDs used as controls"), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Everything done!"), file=logFile, append=TRUE, sep = "\n")

