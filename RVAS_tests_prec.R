######################################################################################
# CRG 
# Hana SUSAK
# date: 01/02/2016
#------------------------------------------------------------------------------------
# SKAT-O, MiST and KBAC test, with or without permutations for controls
#------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calling the script: Rscript germline_risk_test_v2.R --h
#------------------------------------------------------------------------------------
######################################################################################

rm(list=ls())
time.f <- Sys.time()
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MiST))
suppressPackageStartupMessages(library(SKAT))
suppressPackageStartupMessages(library(KBAC,lib.loc="/software/xe/el6.3/R_libs"))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(INLA, lib.loc='/software/xe/el6.3/R_libs/INLA_updated'))

parser <- ArgumentParser()

parser$add_argument("-m", "--mutation.file", type="character", help="input snp/indel file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-d", "--sample_desc_file", type="character", help="samples pca and description file (first column sample ID)", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-f", "--allele_freq", type="character", help="Allele frequency threshold", metavar="value", nargs=1,  default="0.01")
parser$add_argument("-t", "--genes_list_file", type="character", help="Genes (or other analogous list) to test file ", metavar="file", nargs=1)
parser$add_argument("-o", "--out_folder", type="character", help="Folder to set as working directory", metavar="directory", nargs=1)
args <- parser$parse_args()

read.from <- as.character(Sys.getenv("READ_FROM", 'stdin')) 
# read.from <- ''
cat('Standard input is:',read.from,'\n')
# reading/seting global variables
if (read.from=='bash'){
   mc <- as.integer(Sys.getenv("R_MAX_MC_CORES", 1L))    
   seed <- as.integer(Sys.getenv("R_SEED", 160185L)) 
   set.seed(seed)
   INLA_log <- as.logical(Sys.getenv("R_INLA", 'T'))
   MIST_log <- as.logical(Sys.getenv("R_MIST", 'T')) 
   KBAC_log <- as.logical(Sys.getenv("R_KBAC", 'T')) 
   SKATO_log <- as.logical(Sys.getenv("R_SKATO", 'T')) 
   BURDEN_log <- as.logical(Sys.getenv("R_BURDEN", 'T')) 
   maf.log <- as.logical(Sys.getenv("MAF_LOGICAL", 'F')) 
   h.step <- as.numeric(Sys.getenv("H_STEP", 1e-3)) 
} else if (read.from == 'stdin'){ 
   cat("\nEnter a number of cores available: ")
   mc <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nEnter a seed for random slecting of control samples: ")
   seed <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   set.seed(seed)
   cat("\nDo you want to run INLA (T for yes, F for no): ")
   INLA_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nDo you want to run MiST (T for yes, F for no): ")
   if (INLA_log ){
      cat("\nDo you want to use maf as variable for mutations, transformed wiht beta dist.(T for yes, F for no): ")
      maf.log <-  scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
      cat("\nWhat h step for inla you want ot use ? (default 0.001)")
      h.step <- scan(read.from,what=numeric(),nmax=1,quiet=TRUE, nlines=1)
      if(length(h.step)==0){
         h.step<-1e-3
      }
   }
   MIST_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nDo you want to run KBAC (T for yes, F for no): ")
   KBAC_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1)
   cat("\nDo you want to run SKAT-O (T for yes, F for no): ")
   SKATO_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1) 
   cat("\nDo you want to run BURDEN (T for yes, F for no): ")
   BURDEN_log <- scan(read.from,what=logical(),nmax=1,quiet=TRUE, nlines=1) 
} else {
   mc <- as.integer(1L) 
   seed <- as.integer(160185L) 
   set.seed(seed)
   INLA_log <- T
   MIST_log <- T
   KBAC_log <- T
   SKATO_log <- T
   BURDEN_log <- T
   maf.log <- F
   h.step<-1e-3
}

cat("\nNumber of cores useing", mc, " .....\n" )
cat("\nSeed set to", seed, "..... \n" )
cat(paste0("Tests to be performed: ", paste(c('INLA','MiST','KBAC','SKAT-O', 'BURDEN')[c(INLA_log,MIST_log,KBAC_log,SKATO_log, BURDEN_log)], collapse=", ")), "\n")
cat(paste0("INLA maf as variable: ", maf.log, "\n"))
cat(paste0("INLA h step: ", h.step, "\n"))

######################################################################################
# read arguments
######################################################################################

## get command line arguments
mutation.file <- args$mutation.file
samples.info.file <- args$sample_desc_file
af.thr <- as.numeric(args$allele_freq)
genes.file <- args$genes_list_file
out.folder <- args$out_folder
if(FALSE){
   mutation.file <- '/no_backup/GD/projects/RVAS/AR6/mut_data/simulation_1.txt'
   samples.info.file <- '/no_backup/GD/projects/RVAS/AR6/samp_data/samples_info_sim_1.txt'
   af.thr <- 0.01
   
   mutation.file <- '/no_backup/GD/projects/RVAS/QC_CLL_results_folder_2017-03-10-180826/mutations_filtered_samples.txt'
   samples.info.file <- '/no_backup/GD/projects/RVAS/QC_CLL_results_folder_2017-03-10-180826/samples_info_pca.txt'
   af.thr <- 0.01
   genes.file <- '/no_backup/GD/projects/CLL/germ_analysis/results/AF0.005/genes_top_43.txt'
   #genes.file <- NULL
   out.folder <- NULL 
   
   mutation.file <- '/no_backup/GD/projects/RVAS/QC_PCAWG_exome_results_folder_2017-03-08-152715/mutations_filtered_samples.txt'
   samples.info.file <- '/no_backup/GD/projects/RVAS/QC_PCAWG_exome_results_folder_2017-03-08-152715/samples_info_pca.txt'
   af.thr <- 0.01
   genes.file <- '/no_backup/GD/projects/CLL/germ_analysis/data/rahman_genes_fix.txt'
   #genes.file <- NULL
   out.folder <- NULL 
   
   #laura data
   mutation.file <- '/no_backup/GD/projects/HumanDisease/HaplotypeCaller_march15/genotype/genotype_220515/NEW_FEBRUARY_FS/agi35_agi50_nimv3/QC_parseok/QC_results_folder_2017-03-16-144229/mutations_filtered_samples_splicingchange.txt'
   samples.info.file <- '/no_backup/GD/projects/HumanDisease/HaplotypeCaller_march15/genotype/genotype_220515/NEW_FEBRUARY_FS/agi35_agi50_nimv3/QC_parseok/QC_results_folder_2017-03-16-144229/samples_info_pca.txt'
   af.thr <- 0.01
   genes.file <- NULL
   #genes.file <- NULL
   out.folder <- NULL 
   
   mutation.file <- '/no_backup/GD/projects/RVAS/AR1/mut_data/simulation_1.txt'
   samples.info.file <- '/no_backup/GD/projects/RVAS/AR1/samp_data/samples_info_sim_1.txt'
   af.thr <- 0.01
   genes.file <- NULL
   out.folder <- '/users/so/hsusak/' 
   
   mutation.file <- '/users/so/hsusak/Ilnaz_project/QC_results_folder_2017-06-30-173012/mutations_filtered_samples2.txt'
   samples.info.file <- '/users/so/hsusak/Ilnaz_project/QC_results_folder_2017-06-30-173012/samples_info_pca.txt'
   af.thr <- 0.01
   genes.file <- NULL
   out.folder <- '/users/so/hsusak/Ilnaz_project/test' 
   
   
   mutation.file <- '/nfs/no_backup/GD/projects/RVAS/AR1/mut_data/simulation_12.txt'
   samples.info.file <- '/nfs/no_backup/GD/projects/RVAS/AR1/samp_data_t1/samples_info_sim_12/perm_95.txt'
   af.thr <- 0.01
   genes.file <- NULL
   out.folder <- '/users/so/hsusak/' 
   
}

rm(parser)

# normalize paths 
mutation.file <- normalizePath(mutation.file)
samples.info.file <- normalizePath(samples.info.file)
if(!is.null(genes.file)){
   genes.file <- normalizePath(genes.file)   
}
if(!is.null(out.folder)){
   out.folder <- normalizePath(out.folder)   
   cat('Setting to out folder :',out.folder, '\n')
   setwd(out.folder)
} else {
   out.folder <-  paste0('Risk_analysis_',format(time.f, format = "%Y-%m-%d-%H%M%S") )
   cat('Setting to out folder to be :',out.folder, '\n')
   dir.create(out.folder)
   out.folder <- normalizePath(out.folder)   
   setwd(out.folder)
}

cat('Mutations file is:',mutation.file,'\n')
cat('Samples info file is:',samples.info.file,'\n')
cat('Out folder is:',out.folder,'\n')

#creat log file
logFile <-  paste0('Risk_analysis_log_file_',format(time.f, format = "%Y-%m-%d-%H%M%S") , '.txt')
rm(time.f)
cat(timestamp(), file=logFile, append=FALSE, sep = "\n")
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------- Reading input files and variables ------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

cat(paste0("Number of cores: ", mc), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Seed: ", seed), file=logFile, append=TRUE, sep = "\n")

cat(paste0("Tests to be performed: ", paste(c('MiST','KBAC','SKAT-O', 'BURDEN', 'INLA')[c(MIST_log,KBAC_log,SKATO_log, BURDEN_log, INLA_log)], collapse=", ")),
    file=logFile, append=TRUE, sep = "\n")
cat(paste0("INLA maf as variable: ", maf.log),  file=logFile, append=TRUE, sep = "\n")
cat(paste0("INLA h step: ", h.step),  file=logFile, append=TRUE, sep = "\n")

cat(paste0("Mutations file: ", mutation.file), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Samples info file: ", samples.info.file), file=logFile, append=TRUE, sep = "\n")
cat(paste0("AF threshold: ", af.thr), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Genes file: ", genes.file), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Out folder: ", out.folder), file=logFile, append=TRUE, sep = "\n")

######################################################################################
# functions 
######################################################################################

if(INLA_log) {
   doit.inla <- function(y.bin, Xi, G, Z.func, maf, gene){
      M <- cbind(as.matrix(Xi),G%*%as.matrix(Z.func))
      #Qz <- diag(dbeta(maf, weight.beta[1], weight.beta[2])^2)
      if(maf.log){
         weight.beta <- c(1,2)
         Z.func$maf <- diag(dbeta(maf, weight.beta[1], weight.beta[2])^2)   
      }
      
      Y <- y.bin
      n <- nrow(G)
      #M <- M[,-1, drop=F]
      id.z <- 1:n
      formula.null <- Y ~ -1 + as.matrix(Xi)
      Z <- G #%*%as.matrix(Z.func)
      formula <- Y ~ -1 + as.matrix(M) + f(id.z, model="z",  Z=Z) #Cmatrix=Qz,
      result <- inla(formula, family = "logistic", data = list(Y=Y, id.z=1:n, M=M),
                     control.compute=list(dic=T,cpo=TRUE,config=TRUE),num.threads=1,
                     control.inla=list( h=h.step))
      #result2 <- inla(formula, family = "binomial", data = list(Y=Y, id.z=1:n, M=M),control.compute=list(dic=T,cpo=TRUE,config=TRUE),num.threads=1)
      #save(result,file=paste0(gene, '_logistic.RData'))
      #save(result2,file=paste0(gene, '_binomial.RData'))
      
      result.null <- inla(formula.null,family = "logistic",  data = list(Y=Y, id.z=1:n, M=as.matrix(Xi) ),
                          control.compute=list(dic=T,cpo=TRUE,config=TRUE), num.threads=1,
                          control.inla=list( h=h.step))
      #result2.null <- inla(formula.null,family = "binomial",  data = list(Y=Y, id.z=1:n, M=as.matrix(Xi) ),control.compute=list(dic=T,cpo=TRUE,config=TRUE), num.threads=1)
      #save(result.null,file=paste0(gene, '_logistic_null.RData'))
      #save(result2.null,file=paste0(gene, '_binomial_null.RData'))
      
      dic.diff.dic <- result.null$dic$dic - result$dic$dic 
      cpo.diff.cpo <- -mean(log(result.null$cpo$cpo[result.null$cpo$cpo>0]),na.rm=T) - ( -mean(log(result$cpo$cpo[result$cpo$cpo>0]),na.rm=T) ) 
      names(dic.diff.dic) <- 'DIC.diff'
      names(cpo.diff.cpo) <- 'CPO.diff'
      #means <- result$summary.fixed[-c(1:ncol(Xi)),1]
      #names(means) <-  rownames(result$summary.fixed)[-c(1:ncol(Xi))]
      
      #rand.mean <- result$summary.hyperpar[1]
      #names(rand.mean) <-  rownames(result$summary.hyperpar)
      test.ex.temp <- c(dic.diff.dic,cpo.diff.cpo)
      
      return(test.ex.temp)
   } 
   
   test.inla.logit <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, numeric.var, dummy.var, agg.col, pi){   
      write(paste('starting to test gene: ',gene), file=paste('inla_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( agg.col == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         #samp.start <- which(colnames(mutations) == '#LocalAF')+1
         samp.cols <-  which(!grepl('#',colnames(mutations)))
         df <- mutations[,c(gene.col,snpID.col,samp.cols)]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         G <- na.omit(G)
         G <- G[,colSums(G, na.rm=T)!=0 , drop=FALSE]
         if(nrow(G)>0 && ncol(G)>0){          
            cases.samp.nona <- which(rownames(G) %in% cases.samp)
            cases.samp.nona <- rownames(G)[cases.samp.nona]
            controls.samp.nona <- which(rownames(G) %in% controls.samp)
            controls.samp.nona <- rownames(G)[controls.samp.nona]
            tot.nona.cases <- length(cases.samp.nona)
            tot.nona.controls <- length(controls.samp.nona)
            
            mut.total <- ncol(G)
            mut.cases <- sum(colSums(as.matrix(G[cases.samp.nona , ,drop=F]), na.rm=T) !=0)
            mut.controls <- sum(colSums(as.matrix(G[controls.samp.nona, ,drop=F]), na.rm=T) !=0)
            num.cases <- sum(rowSums(as.matrix(G[cases.samp.nona , ,drop=F]), na.rm=T) !=0)
            num.controls <- sum(rowSums(as.matrix(G[controls.samp.nona, ,drop=F]), na.rm=T) !=0)
            
            maf.df <- pop.mut[ colnames(G),]         
            maf <- maf.df$`#LocalAF` + 0.00001
            names(maf)  <- maf.df$`#snpID`
            rm(maf.df)
            
            y.bin <- as.numeric(rownames(G) %in% cases.samp.nona)
            
            if (numeric.var != ""){
               Z.func <- data.frame(intercept_Z=1, temp=mutations[ colnames(G), c(numeric.var)])
               Z.func[is.na(Z.func$temp), 'temp'] <- median(Z.func[, 'temp'], na.rm=T) 
               colnames(Z.func)[-1] <- numeric.var
            } else {
               Z.func <- data.frame(intercept_Z=rep(1, times=ncol(G)))
            }         
            rownames(Z.func) <- colnames(G)
            
            if (dummy.var != ""){
               Z.func2 <- data.frame(temp=as.character(mutations[ colnames(G), c(dummy.var)]))
               for(level in unique(Z.func2$temp)){
                  Z.func2[paste("dummy", level, sep = "_")] <- ifelse(Z.func2$temp == level, 1, 0)
               }        
               rownames(Z.func2) <- colnames(G)
               Z.func2 <- Z.func2[, -c(1:2), drop=F]  
               Z.func <- cbind.data.frame(Z.func[colnames(G),, drop=F], Z.func2[colnames(G),, drop=F] )
               rm(Z.func2)
            } 
            
            Xi <- Xi[rownames(G),, drop=F]
            
            test.ex <- (try( doit.inla(y.bin, Xi, G, Z.func, maf, gene), TRUE))
            #print(test.ex)
            dic.diff <-  as.numeric(test.ex[1])
            cpo.diff <-  as.numeric(test.ex[length(test.ex)])
         } else {
            mut.total <- 0
            mut.cases <- 0
            mut.controls <- 0
            num.cases <-0
            num.controls <- 0
            tot.nona.cases <- 0
            tot.nona.controls <- 0
            dic.diff <-  NA
            cpo.diff <- NA
         }
         
      } else {
         mut.total <- 0
         mut.cases <- 0
         mut.controls <- 0
         num.cases <-0
         num.controls <- 0
         tot.nona.cases <- length(cases.samp)
         tot.nona.controls <-  length(controls.samp)
         dic.diff <-  NA
         cpo.diff <- NA
         
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls,                 
                 no.na.cases=tot.nona.cases, no.na.controls=tot.nona.controls,
                 dic.diff=dic.diff, cpo.diff=cpo.diff ))     ) 
   }
   
   perm.mc.inla <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per, col.snps.numeric, col.snps.dummy, agg.col, pi) {
      write('starting for this permutation!', file=paste('inla_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.inla.logit(x, mutations.per[mutations.per[,agg.col]==x,],
                                                          cases.samples,controls.samples.per,
                                                          pop.mut.df.per, X.per, col.snps.numeric, col.snps.dummy, agg.col,pi))
      
      write('Done with this permutation!', file=paste('inla_perm',pi, sep='_'), sep = "\t", append=T)
      
      df.temp <- as.data.frame(t(df))
      
      #df.temp <- data.frame(df.temp)
      
      colnames(df.temp) <- paste(colnames(df.temp), '_perm',pi, sep='')
      
      write('going out of perm.mc.inla function!', file=paste('inla_perm',pi, sep='_'), sep = "\t", append=T)
      
      return(df.temp)
      
   }
   
}

####

if(MIST_log) {
   doit.mist <- function(y.bin, X, G, Z.func, maf){
      #print(paste('Gene ',gene, 'is tasted'))
      test.ex <- logit.weight.test(y.bin,X,G,Z.func,maf)
      return(test.ex)
   } 
   
   test.mist.logit <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, numeric.var, dummy.var, agg.col, pi){   
      write(paste('starting to test gene: ',gene), file=paste('mist_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( agg.col == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#LocalAF')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         
         G <- t(as.matrix( df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         G <- na.omit(G)
         G <- G[,colSums(G, na.rm=T)!=0 , drop=FALSE]
         if(nrow(G)>0 && ncol(G)>0){          
            cases.samp.nona <- which(rownames(G) %in% cases.samp)
            cases.samp.nona <- rownames(G)[cases.samp.nona]
            controls.samp.nona <- which(rownames(G) %in% controls.samp)
            controls.samp.nona <- rownames(G)[controls.samp.nona]
            tot.nona.cases <- length(cases.samp.nona)
            tot.nona.controls <- length(controls.samp.nona)
            
            mut.total <- ncol(G)
            mut.cases <- sum(colSums(as.matrix(G[cases.samp.nona , ,drop=F]), na.rm=T) !=0)
            mut.controls <- sum(colSums(as.matrix(G[controls.samp.nona, ,drop=F]), na.rm=T) !=0)
            num.cases <- sum(rowSums(as.matrix(G[cases.samp.nona , ,drop=F]), na.rm=T) !=0)
            num.controls <- sum(rowSums(as.matrix(G[controls.samp.nona, ,drop=F]), na.rm=T) !=0)
            
            maf.df <- pop.mut[ colnames(G),]         
            maf <- maf.df$`#LocalAF` + 0.00001
            names(maf)  <- maf.df$`#snpID`
            rm(maf.df)
            
            y.bin <- as.numeric(rownames(G) %in% cases.samp.nona)
            
            if (numeric.var != ""){
               Z.func <- data.frame(intercept=1, temp=mutations[ colnames(G), c(numeric.var)])
               Z.func[is.na(Z.func$temp), 'temp'] <- mean(Z.func[, 'temp'], na.rm=T) 
               colnames(Z.func)[-1] <- numeric.var
            } else {
               Z.func <- data.frame(intercept=rep(1, times=ncol(G)))
            }
            
            rownames(Z.func) <- colnames(G)
            
            if (dummy.var != ""){
               Z.func2 <- data.frame(temp=as.character(mutations[ colnames(G), c(dummy.var)]))
               for(level in unique(Z.func2$temp)){
                  Z.func2[paste("dummy", level, sep = "_")] <- ifelse(Z.func2$temp == level, 1, 0)
               }        
               rownames(Z.func2) <- colnames(G)
               Z.func2 <- Z.func2[, -c(1:2), drop=F]  
               Z.func <- cbind(Z.func[colnames(G),], Z.func2[colnames(G),] )
               rm(Z.func2)
            } 
            
            Xi <- Xi[rownames(G),]
            
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
            tot.nona.cases <- 0
            tot.nona.controls <- 0
            
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
         tot.nona.cases <- length(cases.samp)
         tot.nona.controls <-  length(controls.samp)
         
         p.val.add <-  NA
         p.val.tau <- NA
         p.val.pi <- NA
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls, 
                 no.na.cases=tot.nona.cases, no.na.controls=tot.nona.controls,
                 p.val.pi=p.val.pi,  p.val.tau=p.val.tau,
                 p.val.overall=p.val.add))     ) 
   }
   
   perm.mc.mist <- function(genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per, numeric.var, dummy.var, agg.col, pi) {
      
      write('starting for this permutation!', file=paste('mist_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.mist.logit(gene=x, mutations=mutations.per[mutations.per[,agg.col]==x,],
                                                          cases.samples,controls.samples.per,
                                                          pop.mut.df.per, X.per, numeric.var, dummy.var, agg.col, pi))
      
      #write('Done with this permutation!', file=paste('mist_perm',pi, sep='_'), sep = "\t", append=T)
      
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
   
   test.skat.binary <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, agg.col, pi){   
      write(paste('starting to test gene: ',gene), file=paste('skat0_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( agg.col == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#LocalAF')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         
         mut.total <- ncol(G)
         mut.cases <- sum(colSums(as.matrix(G[cases.samp , ,drop=F]), na.rm=T) !=0)
         mut.controls <- sum(colSums(as.matrix(G[controls.samp, ,drop=F]), na.rm=T) !=0)
         num.cases <- sum(rowSums(as.matrix(G[cases.samp , ,drop=F]), na.rm=T) !=0)
         num.controls <- sum(rowSums(as.matrix(G[controls.samp, ,drop=F]), na.rm=T) !=0)
         tot.nona.cases <- sum(rowSums(is.na(G[cases.samp, ,drop=F]))==0)
         tot.nona.controls <- sum(rowSums(is.na(G[controls.samp, ,drop=F]))==0)
         
         
         maf.df <- pop.mut[ colnames(G),]         
         maf <- maf.df$`#LocalAF` + 0.00001
         names(maf)  <- maf.df$`#snpID`
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
         tot.nona.cases <- length(cases.samp)
         tot.nona.controls <-  length(controls.samp)
         
         p.val <-  NA
         
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls, 
                 no.na.cases=tot.nona.cases, no.na.controls=tot.nona.controls,
                 p.val.overall=p.val))     ) 
   }
   
   perm.mc.skato <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per,agg.col, pi) {
      write('starting for this permutation!', file=paste('skat0_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.skat.binary(x, mutations.per[mutations.per[,agg.col]==x,],
                                                           cases.samples,controls.samples.per,
                                                           pop.mut.df.per, X.per,agg.col,  pi))
      
      write('Done with this permutation!', file=paste('skat0_perm',pi, sep='_'), sep = "\t", append=T)
      
      df.temp <- as.data.frame(t(df))
      
      df.temp <- data.frame(df.temp, p.adjust = p.adjust(df.temp$p.val.overall, "BH", n=length(genes.all)))
      
      colnames(df.temp) <- paste(colnames(df.temp), '_perm',pi, sep='')
      
      write('going out of perm.mc.skato function!', file=paste('skat0_perm',pi, sep='_'), sep = "\t", append=T)
      
      return(df.temp)
      
   }
   
}

####

if(BURDEN_log) {
   doit.burden <- function(y.bin, X, G, maf){
      
      obj.s <- SKAT_Null_Model(y.bin ~ X, out_type="D",Adjustment=FALSE)
      weights <- Get_Logistic_Weights_MAF(maf)   
      test.ex <- SKATBinary(G,obj.s,kernel="linear.weighted",weights=weights,method="Burden")
      
      return(test.ex)
   } 
   
   test.burden.binary <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, agg.col, pi){   
      write(paste('starting to test gene: ',gene), file=paste('burden_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( agg.col == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#LocalAF')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         
         mut.total <- ncol(G)
         mut.cases <- sum(colSums(as.matrix(G[cases.samp , ,drop=F]), na.rm=T) !=0)
         mut.controls <- sum(colSums(as.matrix(G[controls.samp, ,drop=F]), na.rm=T) !=0)
         num.cases <- sum(rowSums(as.matrix(G[cases.samp , ,drop=F]), na.rm=T) !=0)
         num.controls <- sum(rowSums(as.matrix(G[controls.samp, ,drop=F]), na.rm=T) !=0)
         tot.nona.cases <- sum(rowSums(is.na(G[cases.samp, ,drop=F]))==0)
         tot.nona.controls <- sum(rowSums(is.na(G[controls.samp, ,drop=F]))==0)
         
         
         
         maf.df <- pop.mut[ colnames(G),]         
         maf <- maf.df$`#LocalAF` + 0.00001
         names(maf)  <- maf.df$`#snpID`
         rm(maf.df)
         
         y.bin <- as.numeric(rownames(G) %in% cases.samp)
         
         
         #if( (mut.cases>0 & mut.controls>0) | (mut.cases==0 & num.controls>2) | (mut.controls==0 & num.cases>2)){
         test.ex <- (try( doit.burden(y.bin, as.matrix(Xi), G, maf), TRUE))
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
         tot.nona.cases <- length(cases.samp)
         tot.nona.controls <-  length(controls.samp)
         
         p.val <-  NA
         
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,                 
                 num.cases=num.cases, num.controls=num.controls, 
                 no.na.cases=tot.nona.cases, no.na.controls=tot.nona.controls,
                 p.val.overall=p.val))     ) 
   }
   
   perm.mc.burden <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per,agg.col, pi) {
      write('starting for this permutation!', file=paste('burden_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.burden.binary(x, mutations.per[mutations.per[,agg.col]==x,],
                                                             cases.samples,controls.samples.per,
                                                             pop.mut.df.per, X.per,agg.col,  pi))
      
      write('Done with this permutation!', file=paste('burden_perm',pi, sep='_'), sep = "\t", append=T)
      
      df.temp <- as.data.frame(t(df))
      
      df.temp <- data.frame(df.temp, p.adjust = p.adjust(df.temp$p.val.overall, "BH", n=length(genes.all)))
      
      colnames(df.temp) <- paste(colnames(df.temp), '_perm',pi, sep='')
      
      write('going out of perm.mc.burden function!', file=paste('burden_perm',pi, sep='_'), sep = "\t", append=T)
      
      return(df.temp)
      
   }
   
}

####

if(KBAC_log) {
   doit.kbac <- function(y.bin,G ){
      
      casectrl.dat <- as.matrix(cbind(y.bin,G))
      
      alpha <-  2.5e-06 #0.00001
      num.permutation <- 1000000
      quiet <- T
      alternative <- 1
      maf.upper.bound <- 1
      
      test.ex <- KbacTest(casectrl.dat, alpha, num.permutation, quiet, maf.upper.bound, alternative)
      
      return(test.ex)
   } 
   
   test.kbac <- function(gene, mutations, cases.samp, controls.samp, pop.mut, Xi, agg.col,  pi){   
      write(paste('starting to test gene: ',gene), file=paste('kbac_perm',pi, sep='_'), sep = " ", append=T)
      
      if(nrow(mutations) > 0){
         gene.col <- which( agg.col == colnames(mutations))
         snpID.col <- which( '#snpID' == colnames(mutations))
         samp.start <- which(colnames(mutations) == '#LocalAF')+1
         df <- mutations[,c(gene.col,snpID.col,samp.start:ncol(mutations))]
         
         G <- t(as.matrix(df[,-c(1:2)]))
         #colnames(G) <-  df$snpID
         G <- na.omit(G)
         G <- G[,colSums(G, na.rm=T)!=0 , drop=FALSE]
         if(nrow(G)>0 && ncol(G)>0){  
            cases.samp.nona <- which(rownames(G) %in% cases.samp)
            cases.samp.nona <- rownames(G)[cases.samp.nona]
            controls.samp.nona <- which(rownames(G) %in% controls.samp)
            controls.samp.nona <- rownames(G)[controls.samp.nona]
            tot.nona.cases <- length(cases.samp.nona)
            tot.nona.controls <- length(controls.samp.nona)
            
            mut.total <- ncol(G)
            mut.cases <- sum(colSums(as.matrix(G[cases.samp.nona ,  ,drop=F]), na.rm=T) !=0)
            mut.controls <- sum(colSums(as.matrix(G[controls.samp.nona, ,drop=F]), na.rm=T) !=0)
            num.cases <- sum(rowSums(as.matrix(G[cases.samp.nona ,  ,drop=F]), na.rm=T) !=0)
            num.controls <- sum(rowSums(as.matrix(G[controls.samp.nona,  ,drop=F]), na.rm=T) !=0)
            
            maf.df <- pop.mut[ colnames(G),]         
            maf <- maf.df$`#LocalAF` + 0.00001
            names(maf)  <- maf.df$`#snpID`
            rm(maf.df)
            
            y.bin <- as.numeric(rownames(G) %in% cases.samp.nona)
            
            
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
            tot.nona.cases <- 0
            tot.nona.controls <- 0
            p.val <-  NA
         }
         
      } else {
         mut.total <- 0
         mut.cases <- 0
         mut.controls <- 0
         num.cases <-0
         num.controls <- 0
         tot.nona.cases <- length(cases.samp)
         tot.nona.controls <-  length(controls.samp)
         
         p.val <-  NA
      }
      
      
      return((c( total.mut=mut.total , cases.mut=mut.cases, controls.mut=mut.controls,
                 num.cases=num.cases, num.controls=num.controls,
                 no.na.cases=tot.nona.cases, no.na.controls=tot.nona.controls,
                 p.val.overall=p.val))) 
   }
   
   perm.mc.kbac <- function (genes.all, mutations.per, cases.samples, controls.samples.per, pop.mut.df.per, X.per,agg.col, pi) {
      write('starting for this permutation!', file=paste('kbac_perm',pi, sep='_'), sep = "\t", append=T)
      
      df <- sapply(genes.all, function(x) test.kbac(x, mutations.per[mutations.per[,agg.col]==x,],
                                                    cases.samples,controls.samples.per,
                                                    pop.mut.df.per, X.per, agg.col, pi))
      
      write('Done with this permutation!', file=paste('kbac_perm',pi, sep='_'), sep = "\t", append=T)
      
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
#sample.start <- which('#snpID'==colnames(mutations.df))+1

mutations.df <- mutations.df[rowSums(mutations.df[,samples.cols], na.rm=T)!=0 ,]
cat("\nIn input file there are", nrow(mutations.df), "unique snps/indels  \n ....." )

# annotation and samples columnts 
annot.cols <- which(grepl('#',colnames(mutations.df)))
samples.cols <- which(!grepl('#',colnames(mutations.df)))

# samples info file
samples.info <- read.table(file=samples.info.file, header=T, sep="\t", quote="")
cat("\nThere is info for", nrow(samples.info), "samples. \n First column should be unique sample ID, as in multi-sample file \n ....." )

#colnames(samples.info)[1] <- 'samples.my'
samples.ids.col <-  colnames(samples.info)[1] #samples info file first coumn = sample ids

all.samples <- unique(as.character(colnames(mutations.df[,samples.cols])))
if ( all(sapply(all.samples, function(x) x %in% samples.info[,samples.ids.col])) ){
   samples.info <- samples.info[samples.info[,samples.ids.col] %in% all.samples,]  
} else {
   stop("There must be info for all samples from mutation file. Some samples are mismached or missing.")
}

########################
# X matrix infos
########################

if (read.from=='bash'){
   pca.cols <- Sys.getenv("R_COV_MAT_COLS", unset="")  
   if(pca.cols != ""){
      pca.cols <- as.character(pca.cols)
      pca.cols <- trimws(unlist(strsplit(pca.cols, ',')))
      #pca.cols <- as.numeric(pca.cols)
      x <- which( colnames(samples.info) %in% pca.cols)
      if(length(x)==length(pca.cols)){
         pca.x <- samples.info[,c(samples.ids.col,pca.cols)]
         rownames(pca.x) <- pca.x[,samples.ids.col]
         if(any(sapply(pca.x[-1], class) != 'numeric')){
            stop(paste("One of the columns for covariates (",paste(pca.cols, collapse=', ' ),') is not numeric column'))            
         }            
      } else {
         stop(paste("One of the columns for covariates (",paste(pca.cols, collapse=', ' ),
                    ') is not present in file',samples.info.file))
      }      
      rm(x)
   } else {
      cat("\nAs R_COV_MAT_COLS is not set, covariate matrix will not be used in further analysis")
   }
   rm(pca.cols)  
} else if (read.from == 'stdin'){ 
   cat("\n ---------------------------------------------------------------------------------------------\n")
   cat("\n ------------------------------------ Chooseing covariates -----------------------------------\n")
   cat("\n ---------------------------------------------------------------------------------------------\n")
   
   cols.types <- sapply(samples.info, class)
   num.cols <- cols.types[cols.types =='numeric' |cols.types =='integer' ]
   #cat(length(num.cols))
   if(length(num.cols) >0) {
      num.cols <- names(num.cols)
      rm(cols.types)
      
      cat("\nThere are numeric columns: ", 
          paste('\n',1:length(num.cols),'.',num.cols, sep=''),
          "\nWhich one you want to keep included in X (cov mat) matrix (if none, press ENTER): ")
      x <- scan(read.from,what=integer(),nmax=length(num.cols),quiet=TRUE, nlines=1)
      if(length(x)>0){
         x <- num.cols[x]
         pca.x <- samples.info[,c(samples.ids.col,x)]
         rownames(pca.x) <- pca.x[,samples.ids.col]
      }
      rm(x)  
   }   
   rm(num.cols)
} 

if(exists('pca.x')){
   pca.x <- pca.x[pca.x[,samples.ids.col] %in% all.samples , ]  
   cat(paste0("Columns used for covariates are: ", paste0(colnames(pca.x)[-1], collapse=', ')), file=logFile, append=TRUE, sep = "\n")
   cat(paste0("Columns used for covariates are: ", paste0(colnames(pca.x)[-1], collapse=', ')), "\n")
} else {
   cat(paste0("Covariates will not be used"), file=logFile, append=TRUE, sep = "\n")
   cat("Covariates will not be used\n")
   pos.n <- which(samples.ids.col == colnames(samples.info))
   pca.x <- as.data.frame(samples.info[,samples.ids.col, drop=FALSE])
   colnames(pca.x) <- colnames(samples.info)[pos.n]
   #pca.x$dummy <- 1
   rownames(pca.x) <- pca.x[,samples.ids.col]
   rm(pos.n)
}

########################
# AF threshold 
print(paste('\nAllele frequency threshold is: ',af.thr))  

########################
# choose aggregate column
########################
if (read.from=='bash'){
   agg.col <- as.character(Sys.getenv("R_AGG_COL"))  
   if(agg.col==""){
      stop(paste('You need to provide column for aggregation (R_AGG_COL =',agg.col,') from file', mutation.file))
   }
} else if (read.from == 'stdin'){
   cols.types <- sapply(mutations.df[,annot.cols], class)
   car.cols <- cols.types[cols.types =='character' | cols.types =='factor' ]
   car.cols <- names(car.cols)
   rm(cols.types)
   
   cat("\nThere are following character columns: ", 
       paste('\n',1:length(car.cols),'.',car.cols, sep=''),
       "\nWhich one you want to use as aggregation column (e.g. Genes): ")
   agg.col <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   if(length(agg.col)>0){
      agg.col <- car.cols[agg.col]
   } else {
      stop(paste('There is no column for aggregation in file', mutation.file))
   }
   rm(car.cols)
} 

if (!agg.col %in% colnames(mutations.df)[annot.cols]){
   stop('Column you choose to aggregate on has to exist in mutations file!')
}

cat(paste0("Aggregation column: ",agg.col), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Aggregation column: ",agg.col), "\n")

########################
# read genes to test
########################

if(!is.null(genes.file)){
   print(paste('Genes (or other aggregate unit) list file name: ',genes.file))  
   genes.all <- read.table(file=genes.file, header=F, comment.char="", quote="", sep='\t', stringsAsFactors=F, check.names=F)  
   genes.all <- genes.all[!genes.all %in% c('Gene','gene','Gene_Symbol', '#Gene')] # in case of a header
   if (is.data.frame(genes.all)) {
      genes.all <- unique(as.character(genes.all[,1]))
   }
   mutations.df <- mutations.df[mutations.df[,agg.col] %in% genes.all,]
} else {
   genes.all <- unique(as.character(mutations.df[,agg.col]))  
}   
genes.all <- as.character(sort(genes.all))

cat("\nWe will test for : ", length(genes.all), "genes (or aggregetion units)\n ...")
cat("\nWe will use in total: ", nrow(pca.x), "patients (some cases, some controls, some as AF estimators )\n ...")
cat(paste0("Number of Genes to be tested: ", length(genes.all)), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Number of Samples to be tested (case+control): ",  nrow(pca.x)), file=logFile, append=TRUE, sep = "\n")

######################################################################################
# Samples for calculating local population AF
######################################################################################

###############################
# read cases + contorles column
###############################

if (read.from=='bash'){
   proj.col <- as.character(Sys.getenv("R_PROJECT_COL"))    
   cases.name <- as.character(Sys.getenv("R_CASES_NAME"))  
   cases.name <- trimws(unlist(strsplit(cases.name, ',')))
   controls.name <- as.character(Sys.getenv("R_CONTROLS_NAME", unset=""))
   controls.name <- trimws(unlist(strsplit(controls.name, ',')))
   x <- which(proj.col == colnames(samples.info))
   if(length(x)==1){
      options <- as.character(unique(samples.info[,x]))
      cases.names.num <- which(options %in% cases.name )
      if(length(cases.names.num) ==0){
         stop('R_CASES_NAME variable have to match at least one value in R_PROJECT_COL column ')        
      }   
      if(length(controls.name)>0){
         control.names.num <- which(options %in% controls.name )
         if(length(control.names.num) ==0){
            #if there is no matchin controls defined by user, take everything different then cases
            controls.name <- options[!options %in% cases.name]    
         } 
      } else {
         controls.name <- options[!options %in% cases.name]
      }
      rm(options)
   } else {
      stop('R_PROJECT_COL variable has to match excatly one of the character columns names in samples info file')
   }
   rm(x)
} else if (read.from == 'stdin'){ 
   cols.types <- sapply(samples.info, class)
   car.cols <- cols.types[cols.types =='character' | cols.types =='factor' ]
   car.cols <- names(car.cols)
   rm(cols.types)
   
   cat("\n\nWhich column represents cases/controls (Projects)):", 
       paste('\n',1:length(car.cols),'.',car.cols, sep=''))
   cat("\nEnter column number: ")
   x <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   x <- x[x<= length(car.cols)]
   if (length(x)>0){
      proj.col <- car.cols[x]
   } else {
      stop('You have to choose valid column to separete cases and controls')
   }
   rm(x)
   cat("\nColumn(s) you chose is:",proj.col, "\n.....")
   options <- as.character(unique(samples.info[,proj.col]))
   options[is.na(options)] <- 'NA_Missing_info'
   cat("\n\nIn column",proj.col,"there are folowing fields:")
   counts.tb <- table(samples.info[,proj.col], useNA='always')
   names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
   counts.tb <- counts.tb[options]
   cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' samples)','\n', sep=''))
   #cases
   cat("\nEnter  numbers which are considered as cases: ")
   cases.names.num <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1) 
   cases.names.num <- unique(cases.names.num)
   cases.names.num <- cases.names.num[cases.names.num <= length(options)]
   if(length(cases.names.num)>0){   
      cases.name <- options[cases.names.num]
   } else {
      stop('You have to choose valid option which represent cases')    
   }   
   cat("\nSamples whos", proj.col,  "is/are:", paste(cases.name, collapse=', '),
       "are considered as cases. \n.....") 
   #controls
   cat("\nEnter  numbers which are considered as controls (or enter for all others as controls): ")
   controls.names.num <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1) 
   controls.names.num <- unique(controls.names.num)
   controls.names.num <- controls.names.num[controls.names.num <= length(options)]
   controls.names.num <- controls.names.num[! controls.names.num %in% cases.names.num]
   if(length(controls.names.num)>0){   
      controls.name <- options[controls.names.num]
   } else {
      controls.name <- options[! options %in% cases.name]   
   }
   cat("\nSamples whos", proj.col,  "is/are:", paste(controls.name, collapse=', '),
       "are considered as controls. \n.....") 
   
   rm(cases.names.num)
   rm(controls.names.num)
   rm(counts.tb)
   rm(options)
   rm(car.cols)
} 

cat(paste0("Cases defined as: ", paste(cases.name, collapse=', '), " in column ",proj.col,
           " from file ", samples.info.file), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Cases defined as: ", paste(controls.name, collapse=', '), " in column ",proj.col,
           " from file ", samples.info.file), file=logFile, append=TRUE, sep = "\n")


cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('---------------------------------------- Preparation -----------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

# - split cases, controls and AF-set
cases.samples <- as.character(samples.info[samples.info[,proj.col] %in% cases.name, samples.ids.col ])
controls.samples.all <- as.character(samples.info[samples.info[,proj.col] %in% controls.name, samples.ids.col ])

cat(paste("There are in total", length(cases.samples), 'cases'), file=logFile, append=TRUE, sep = "\n")
cat(paste("There are in total", length(controls.samples.all), 'controls'), file=logFile, append=TRUE, sep = "\n")

#clean from popsible unused samples
annot.cols <- which(grepl('#',colnames(mutations.df)))
samples.cases.cols <- which(colnames(mutations.df) %in% cases.samples)
samples.controls.cols <- which(colnames(mutations.df) %in% controls.samples.all)
mutations.df <- mutations.df[,c(annot.cols, samples.cases.cols, samples.controls.cols)]
annot.cols <- which(grepl('#',colnames(mutations.df)))
samples.cols <- which(!grepl('#',colnames(mutations.df)))
rm(samples.cases.cols)
rm(samples.controls.cols)
mutations.df <- mutations.df[rowSums(mutations.df[,samples.cols], na.rm=T)!=0 ,]

# chooce number of permutations for selecting cases and AF samples
if (read.from=='bash'){
   permut <- as.integer(Sys.getenv("R_PERMUTATIONS", 100L))    
   if ( permut < 1){
      stop("You should have number of permutations larger then 0")  
   }
   
   # character columns annotation filtration
   filt.char.col <- as.character(Sys.getenv("R_CHAR_FILT_COL", unset=""))      
   filt.char.val <- as.character(Sys.getenv("R_CHAR_FILT_VAL", unset=""))    
   filt.char <- !((filt.char.col=="") | (filt.char.val == ""))
   if(filt.char){
      filt.char.col <- trimws(unlist(strsplit(filt.char.col, '\\|')))
      filt.char.val <- trimws(unlist(strsplit(filt.char.val, '\\|')))
      filt.char.val <- sapply(filt.char.val, function(x) trimws(unlist(strsplit(x, ','))))
      if(length(filt.char.col)==length(filt.char.val)){
         names(filt.char.val) <- filt.char.col
      } else {
         #log file
         cat(paste("Something was not speciffied correctly and mutations based on any character annotation column - program stoped"),  
             file=logFile, append=TRUE, sep = "\n")        
         stop("Something was not speciffied correctly when filtering based on any character annotation column!")
      }
      for(char.temp in filt.char.col){
         char.val.temp <- filt.char.val[[char.temp]]
         if ((char.temp %in% colnames(mutations.df)[annot.cols]) && 
                (is.character(mutations.df[,char.temp]) || is.factor(mutations.df[,char.temp]))){        
            if( any('NA'== char.val.temp)){
               mutations.df <-  mutations.df[!is.na(as.character(mutations.df[,char.temp])) , ]
            }
            mutations.df <- mutations.df[! mutations.df[,char.temp] %in% char.val.temp | is.na(mutations.df[,char.temp]), ]
            cat("\nMutations whos", char.temp,"is in:", paste(char.val.temp, collapse=', '),"- will be excluded")        
            #log file - snps character annot keept
            options <- as.character(unique(mutations.df[,char.temp]))
            #options[is.na(options)] <- 'NA_Missing_info'
            cat(paste("Mutations are keept if",char.temp,"is matching:",paste(options, collapse=', ' ), sep=' '), 
                file=logFile, append=TRUE, sep = "\n") 
            rm(options)
         } else {
            cat(paste("\nSomething was not speciffied correctly and mutations based on ", char.temp, " annotation column will not be filtred.\n....."))
            #log file
            cat(paste("Something was not speciffied correctly and mutations based on ", char.temp, "  annotation column will not be filtred."),  file=logFile, append=TRUE, sep = "\n")        
         }     
      }
      rm(char.temp)
      rm(char.val.temp)
   } else {
      cat("\n Nothing will be filtered based on character annotations of mutations ")
      cat(paste("Nothing will be filtered based on character annotations of mutations."), file=logFile, append=TRUE, sep = "\n")   
   }
   rm(filt.char)
   rm(filt.char.col)
   rm(filt.char.val)
   
   #numeric columns annotation filtration
   filt.num.col <- as.character(Sys.getenv("R_NUM_FILT_COL", unset=""))    
   filt.num.val <- as.character(Sys.getenv("R_NUM_FILT_VAL", unset=""))  
   filt.num.trh <- as.character(Sys.getenv("R_NUM_FILT_TRH", unset= ""))    
   filt.num <- !((filt.num.col=="") | (filt.num.val=="") | (filt.num.trh==""))
   if(filt.num){
      filt.num.col <- trimws(unlist(strsplit(filt.num.col, '\\|')))
      filt.num.val <- trimws(unlist(strsplit(filt.num.val, '\\|')))
      filt.num.val <- sapply(filt.num.val, function(x) tolower(substr(x,1,1)))
      filt.num.trh <- trimws(unlist(strsplit(filt.num.trh, '\\|')))
      filt.num.trh <- sapply(filt.num.trh, as.numeric)
      
      if(length(filt.num.col)==length(filt.num.val) && length(filt.num.col)==length(filt.num.trh)){
         names(filt.num.val) <- filt.num.col
         names(filt.num.trh) <- filt.num.col
      } else {
         #log file
         cat(paste("Something was not speciffied correctly and mutations based on any numeric annotation column - program stoped"),  
             file=logFile, append=TRUE, sep = "\n")        
         stop("Something was not speciffied correctly when filtering based on any numeric annotation column!")     
      }
      for(num.temp in filt.num.col){
         num.val.temp <- filt.num.val[[num.temp]]
         num.trh.temp <- filt.num.trh[[num.temp]]
         if ((num.temp %in% colnames(mutations.df)[annot.cols]) && is.numeric(mutations.df[,num.temp])){ 
            NA.temp.col <- sum(is.na(mutations.df[,num.temp]))
            if(NA.temp.col >0 ){
               cat(paste(NA.temp.col,"mutations with NA in columns",num.temp,"that are kept.", sep=' '),  file=logFile, append=TRUE, sep = "\n")
            }
            rm(NA.temp.col)
            if(num.val.temp == 'l'){
               mutations.df <- mutations.df[mutations.df[,num.temp] >= num.trh.temp | is.na(mutations.df[,num.temp]) ,  ]
               cat('Everything lower then ',num.trh.temp,' in column ', num.temp, ' will be removed\n')
               #log file - snps numeric annot keept
               cat(paste("Mutations are keept if",num.temp,"is greater or equal then:", num.trh.temp, sep=' '), 
                   file=logFile, append=TRUE, sep = "\n")    
            } else if (num.val.temp == 'h') {
               mutations.df <- mutations.df[mutations.df[,num.temp] <= num.trh.temp | is.na(mutations.df[,num.temp]),  ]
               cat('Everything higher then ',num.trh.temp,' in column ', num.temp, ' will be removed\n')
               #log file - snps numeric annot keept
               cat(paste("Mutations are keept if",num.temp,"is lower or equal then:", num.trh.temp, sep=' '), 
                   file=logFile, append=TRUE, sep = "\n")    
            } else {
               cat('You did not put l or h for column ', num.temp, '. Mutations will not be filtered by this numeric column\n')         
               #log file - snps numeric annot keept
               cat(paste("No mutations is excluded based on",num.temp," - you did not choose correctly, h or l.", sep=' '), 
                   file=logFile, append=TRUE, sep = "\n")    
            }
            
         } else {
            cat("\nSomething was not speciffied correctly and mutations based on any numeric annotation column will not be filtred.\n.....")
            #log file
            cat(paste("Something was not speciffied correctly and mutations based on any numeric annotation column will not be filtred."),  file=logFile, append=TRUE, sep = "\n")         
         }
      }
      rm(num.temp)
      rm(num.val.temp)
      rm(num.trh.temp)
   } else {
      cat("\n Nothing will be filtered based on numeric annotation of mutations ")
      cat(paste("Nothing will be filtered based on numeric annotation of mutations."), file=logFile, append=TRUE, sep = "\n") 
   }
   rm(filt.num.val)
   rm(filt.num.col)
   rm(filt.num.trh)
   rm(filt.num)
} else if (read.from == 'stdin'){ 
   cat("\nHow many controls permutations you want: ")
   permut <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   if ( permut < 1){
      stop("You should enter integer number larger then 0")  
   }  
   cat('\n----------------------------------------------------------------------------------------------\n')
   cat('----------------------------------------- Filtration -----------------------------------------\n')
   cat('----------------------------------------------------------------------------------------------\n') 
   
   ### filter by num and char columns
   cols.types <- sapply(mutations.df[,annot.cols], class)
   char.cols <- cols.types[cols.types =='character' | cols.types =='factor']
   num.cols <- cols.types[cols.types =='numeric' |cols.types =='integer' ]
   rm(cols.types)
   
   #character columns
   cat("\nYou can choose chaarachter  colums to filter by before testin:", 
       paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
   cat("\nEnter the column numbers, space separated: ")
   x <- scan(read.from,what=integer(),nmax=length(char.cols),quiet=TRUE, nlines=1)
   x <- unique(x)
   x <- x[x<= length(char.cols)]       
   if (length(x) <= 0){
      cat("\nYou or choose not to filter or wrong columns \n.....")
      #log file
      cat(paste("You choose not to filter mutations based on any character annotation column."),  file=logFile, append=TRUE, sep = "\n")
   } else {
      cat("\nColumn(s) you choose are:",paste(names(char.cols)[x], collapse=', ' ), "\n.....")
      #log file - snps character annot columns to filter on
      cat(paste("Character column(s) you choose to filter on:",paste(names(char.cols)[x], collapse=', ' )),  file=logFile, append=TRUE, sep = "\n")
      for (i in x){
         temp.col <- names(char.cols[i])
         options <- as.character(unique(mutations.df[,temp.col]))
         options[is.na(options)] <- 'NA_Missing_info'
         cat("\n\nIn column",names(char.cols)[i],"there are folowing fields:")
         counts.tb <- table(mutations.df[,temp.col], useNA='always')
         names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
         counts.tb <- counts.tb[options]
         cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' mutations)','\n', sep=''))
         cat("Enter the numbers which are in front of the options you DON'T want to keep (space separated): ")
         exclude.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
         exclude.options <- unique(exclude.options)
         exclude.options <- exclude.options[exclude.options <= length(options)]
         if (length(exclude.options) > 0) {
            if( 'NA_Missing_info' %in% options[exclude.options]){
               mutations.df <-  mutations.df[!is.na(as.character(mutations.df[,temp.col])) , ]
            }
            mutations.df <- mutations.df[! mutations.df[,temp.col] %in% options[exclude.options]  |  is.na(mutations.df[,temp.col]),  ]
            cat("\nMutations whos", temp.col,"is in:", paste(options[exclude.options], collapse=', '),"- will be excluded")  
            #log file - snps character annot keept
            cat(paste("Mutations are keept if",temp.col,"is matching:",paste(options[-exclude.options], collapse=', ' ), sep=' '), 
                file=logFile, append=TRUE, sep = "\n")            
         }   else {
            cat("\n Nothing from column",temp.col,"will be excluded ")
            cat(paste("Nothing from column",temp.col,"will be excluded.", sep=' '), file=logFile, append=TRUE, sep = "\n") 
         }     
         cat("\n ---------------------------------------------------------------------------------------------")
      }
      rm(exclude.options)         
      rm(options)
      rm(counts.tb)
      rm(i)
      rm(temp.col)
   }  
   rm(x)
   rm(char.cols)
   
   #numeric cols
   cat("\n\nYou can choose numeric colums to filter by:", 
       paste('\n',1:length(num.cols),'.',names(num.cols), sep='' ))
   cat("\nEnter the column numbers, space separated: ")
   x <- scan(read.from,what=integer(),nmax=length(num.cols),quiet=TRUE, nlines=1)
   x <- unique(x)
   x <- x[x<= length(num.cols)]       
   if (length(x) <= 0){
      cat("\nYou  choose not to filter or wrong columns \n.....")
      #log file
      cat(paste("You choose not to filter mutations based on any numeric annotation column."),  file=logFile, append=TRUE, sep = "\n")
   } else {
      cat("\nColumn(s) you choose is:",paste(names(num.cols)[x], collapse=', ' ), "\n.....")
      #log file - snps character annot columns to filter on
      cat(paste("Numeric column(s) you choose to filter on:",paste(names(num.cols)[x], collapse=', ' )),  file=logFile, append=TRUE, sep = "\n")
      for (i in x){
         temp.col <- names(num.cols[i])
         quant.temp.col <- quantile(mutations.df[,temp.col], na.rm=T)
         NA.temp.col <- sum(is.na(mutations.df[,temp.col]))
         if(NA.temp.col >0 ){
            cat("\nIn column",temp.col,"therea are",NA.temp.col,'Na values')
            cat("\n\nDo you want to filter out these mutations (answer y or n)? \nAnswer:")
            NA.ind <- scan(read.from,what=character(),n=1,quiet=TRUE)
            if (NA.ind=='y'){
               mutations.df <- mutations.df[!is.na(mutations.df[,temp.col]),  ]
               cat(paste(NA.temp.col,"mutations with NA in columns",temp.col,"are removed", sep=' '),  file=logFile, append=TRUE, sep = "\n")
            }
            rm(NA.ind)         
            
         }
         cat("\nIn column",temp.col,"quantiles values are:", 
             paste(paste0('\n',quant.temp.col), names(quant.temp.col), sep=' at qunatile ' ) )
         cat("\n\nDo you want to filter out higher or lower then treshold (h for higer, l for lower)? \nAnswer: ")
         num.ind <- scan(read.from,what=character(),n=1,quiet=TRUE)
         num.ind <- tolower(num.ind)
         num.ind <- substr(num.ind,1,1)
         cat("\nWhat is treshold you want to use: ")
         trehsold <- scan(read.from,what=numeric(),n=1,quiet=TRUE)       
         if(num.ind == 'l'){
            mutations.df <- mutations.df[mutations.df[,temp.col] >= trehsold |  is.na(mutations.df[,temp.col]),  ]
            message('Everything lower then ',trehsold,' in column ', temp.col, ' will be removed')
            #log file - snps numeric annot keept
            cat(paste("Mutations are keept if",temp.col,"is greater or equal then:", trehsold, sep=' '), 
                file=logFile, append=TRUE, sep = "\n")    
         } else if (num.ind == 'h') {
            mutations.df <- mutations.df[mutations.df[,temp.col] <= trehsold |  is.na(mutations.df[,temp.col]),  ]
            message('Everything higher then ',trehsold,' in column ', temp.col, ' will be removed')
            #log file - snps numeric annot keept
            cat(paste("Mutations are keept if",temp.col,"is lower or equal then:", trehsold, sep=' '), 
                file=logFile, append=TRUE, sep = "\n")    
         } else {
            message('You did not put l or h. It will not be filtered by this column')         
            #log file - snps numeric annot keept
            cat(paste("No mutations is excluded based on",temp.col," - you did not choose correctly, h or l.", trehsold, sep=' '), 
                file=logFile, append=TRUE, sep = "\n")    
         }
         
         cat("\n ---------------------------------------------------------------------------------------------\n")
      }    
      rm(NA.temp.col)
      rm(quant.temp.col)
      rm(num.ind)
      rm(trehsold)
      rm(i)
      rm(temp.col)   
   }  
   rm(x)   
   rm(num.cols)   
} else {
   permut <- 100
}


#################################################################################
# Collapsing categoric column to dummy variable (for mist and and ina)
#################################################################################

if (MIST_log | INLA_log){
   cat("\n ---------------------------------------------------------------------------------------------\n")
   cat("\n --------------------------- Random part variables, for mutations ----------------------------\n")
   cat("\n --------------------- Numeric and Dummy variable selection and grouping ---------------------\n")
   cat("\n ---------------------------------------------------------------------------------------------\n")
   
   # lof file
   cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
   cat('---------------------        Random part variables, for mutations       ----------------------', file=logFile, append=TRUE, sep = "\n")
   cat('--------------------- Numeric and Dummy variable selection and grouping ----------------------', file=logFile, append=TRUE, sep = "\n")
   cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")
   
   if (read.from=='bash'){
      col.snps.numeric <- as.character(Sys.getenv("R_NUMERIC_VAR", unset=""))   
      
      col.snps.dummy <- as.character(Sys.getenv("R_DUMMY_VAR", unset=""))    
      combine.options <- as.character(Sys.getenv("R_DUMMY_MERGE", unset=""))  
      rename.options <- as.character(Sys.getenv("R_DUMMY_RENAME", unset="")) 
      
      combine.options <- trimws(unlist(strsplit(combine.options, '\\|')))
      combine.options <- lapply(combine.options, function(x) trimws(unlist(strsplit(x, ','))))
      rename.options <- trimws(unlist(strsplit(rename.options, '\\|')))
      cat(paste0("\nColumn you choose to be used as character dummy is:",paste(col.snps.dummy, collapse=', ' )), file=logFile, append=TRUE, sep = "\n")
      cat(paste0("\nColumn you choose to be used as numeric:",paste(col.snps.numeric, collapse=', ' )),file=logFile, append=TRUE, sep = "\n")
      
      if (length(combine.options) > 0 ){
         for (i.d in 1:length(combine.options)){
            comb.opt <- combine.options[[i.d]]
            rename.opt <- rename.options[i.d]
            if( 'NA_Missing_info' %in% comb.opt){
               mutations.df[is.na(mutations.df[,col.snps.dummy]),col.snps.dummy] <- rename.opt   
               if (length(combine.options) > 1){
                  mutations.df[mutations.df[,col.snps.dummy] %in% comb.opt, col.snps.dummy]<- rename.opt
               }
            } else {
               mutations.df[mutations.df[,col.snps.dummy] %in% comb.opt &
                               !is.na(mutations.df[,col.snps.dummy]), col.snps.dummy]<- rename.opt      
            }
            cat(paste("Mutations which have", col.snps.dummy,"matching:",paste0(comb.opt,collapse=', '),
                      " - will be renamed to",rename.opt), 
                file=logFile, append=TRUE, sep = "\n")
         }
         rm(comb.opt)
         rm(rename.opt)
         rm(i.d)
         
      }
      rm(rename.options)
      rm(combine.options)
   } else if (read.from == 'stdin'){ 
      annot.cols <- which(grepl('#',colnames(mutations.df)) )
      samples.cols <- which(!grepl('#',colnames(mutations.df)))
      
      cols.types <- sapply(mutations.df[,annot.cols], class)
      char.cols <- cols.types[cols.types =='character']
      num.cols <- cols.types[cols.types =='numeric' | cols.types =='integer']
      rm(cols.types)
      cat("\n\nYou can choose character column to be used as dummy variable (for random part - snp vise) :", 
          paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
      cat("\nEnter the column number (or if you don't want to do this, press enter): ")
      x.col <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
      x.col <- unique(x.col)
      x.col <- x.col[x.col<= length(char.cols)] 
      
      if (length(x.col) <= 0){
         cat("\nYou choose not to create dummy variable for snps\n.....")
         cat('You choose not to create dummy variable for snps. ', file=logFile, append=TRUE, sep = "\n")
         col.snps.dummy <- ""
      } else {
         col.snps.dummy <- names(char.cols[x.col])
         cat("\nColumn you choose to be used as dummy is:",paste(col.snps.dummy, collapse=', ' ), "\n.....")
         cat(paste0('Column you choose to be used as dummy is: ',col.snps.dummy), file=logFile, append=TRUE, sep = "\n")
         options <- as.character(unique(mutations.df[,col.snps.dummy]))
         options[is.na(options)] <- 'NA_Missing_info'
         cat("\n\nIn column",col.snps.dummy,"there are folowing fields:")
         counts.tb <- table(mutations.df[,col.snps.dummy], useNA='always')
         names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
         counts.tb <- counts.tb[options]
         cat("\n",paste(1:length(options),".", options ," (", counts.tb ," mutations)","\n", sep=''))
         cat("Do you want to merge and  give new name to some of the column's values, for dummy variable? Answer y or n : ")
         while( !exists('col.ind') ||  length(col.ind)==0 || !col.ind %in% c('y','n')  ){
            col.ind <- try( scan(read.from,what=character(),n=1,quiet=TRUE, nlines=1) ,  silent=T )  
            col.ind <- tolower(substr(col.ind,1,1))
         }
         
         while(col.ind=='y'){
            cat("Which one you want to combine: ")
            combine.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
            combine.options <- unique(combine.options)
            combine.options <- combine.options[combine.options <= length(options)]
            cat("How you would like to name them: ")
            rename.options <- scan(read.from,what=character(),nmax=1,quiet=TRUE, nlines=1)
            if( 'NA_Missing_info' %in% options[combine.options]){
               mutations.df[is.na(mutations.df[,col.snps.dummy]),col.snps.dummy] <- rename.options   
               if (length(combine.options) > 1){
                  mutations.df[mutations.df[,col.snps.dummy] %in% options[combine.options], col.snps.dummy]<- rename.options
               }
            } else {
               mutations.df[mutations.df[,col.snps.dummy] %in% options[combine.options] &
                               !is.na(mutations.df[,col.snps.dummy]), col.snps.dummy]<- rename.options      
            }
            cat(paste("Mutations which have", col.snps.dummy,"matching:",paste0(options[combine.options],collapse=', '),
                      " - will be renamed to",rename.options), 
                file=logFile, append=TRUE, sep = "\n")
            
            ##
            options <- as.character(unique(mutations.df[,col.snps.dummy]))
            options[is.na(options)] <- 'NA_Missing_info'
            if(length(options)>1){
               cat("\n\nIn column",col.snps.dummy,"there are folowing fields now:")
               counts.tb <- table(mutations.df[,col.snps.dummy], useNA='always')
               names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
               counts.tb <- counts.tb[options]
               cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' mutations)','\n', sep=''))
               cat("Do you want to merge more? Answer y or n : ")
               rm(col.ind )
               while( !exists('col.ind') || length(col.ind)==0 || !col.ind %in% c('y','n')  ){
                  col.ind <- try( scan(read.from,what=character(),n=1,quiet=TRUE, nlines=1) ,  silent=T)  
                  col.ind <- tolower(substr(col.ind,1,1))
               }               
            } else {
               col.ind <- 'n'
            }
            rm(combine.options)
            rm(rename.options)
         }  
         rm(counts.tb)   
         rm(col.ind)
         rm(options)
      }
      rm(char.cols)
      rm(x.col) 
      
      cat("\n\nYou can choose numeric column to be used as variable (for random part - snp vise) :", 
          paste('\n',1:length(num.cols),'.',names(num.cols), sep='' ))
      cat("\nEnter the column number (or if you don't want to do this, press enter): ")
      x.col <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
      x.col <- unique(x.col)
      x.col <- x.col[x.col<= length(num.cols)] 
      if (length(x.col) <= 0){
         cat("\nYou choose not to create numeric variable for snps \n.....")
         cat('You choose not to create numeric variable for snps. ', file=logFile, append=TRUE, sep = "\n")
         col.snps.numeric <- ""
      } else {
         col.snps.numeric <- names(num.cols[x.col])
      }
      rm(num.cols)
      rm(x.col) 
   } 
}
#rm(char.cols)

######################################################################################
# preparing and looping
######################################################################################

annot.cols <- which(grepl('#',colnames(mutations.df)) )  
# - calculate local AF on all controls samples
pop.cols <- which(colnames(mutations.df) %in% controls.samples.all)
pop.mut <- mutations.df[,c( annot.cols , pop.cols)]
pop.annot.cols <- which(grepl('#',colnames(pop.mut)) )  
rownames(pop.mut) <- pop.mut$`#snpID`
#freq <- rowSums(pop.mut[,-pop.annot.cols], na.rm=T)/(2*length(controls.samples.all))
freq <- rowSums(pop.mut[,-pop.annot.cols], na.rm=T)/(2*apply(pop.mut[,-pop.annot.cols], 1 , function(x) sum(!is.na(x))  ))

pop.mut$`#LocalAF` <- sapply(freq, function(x) min(x, 1-x) ) # corrected minor af        
pop.mut.df.all.contr <- pop.mut[,c('#snpID', '#LocalAF')]
#colnames(pop.mut.df.all.contr) <- c('#snpID', '#LocalAF')
rm(pop.cols)
rm(pop.mut)
rm(freq)    
rm(pop.annot.cols)  
# SAVE this pop.mut.df.all.contr table

#df.all <- data.frame()
same.perm <- FALSE

af.samples <- list()
controls.samples <- list()
pop.mut.df <- list()  

for (p in 1:ceiling(permut/mc)){
   pp <- ((p-1)*mc+1):(p*mc)
   cat("\nPermutaion number:",pp,"!\n.....")
   # split controls- for local AF estimation
   if (p==1 | !same.perm ){
      if (read.from=='bash'){
         controls.num <- as.integer(Sys.getenv("R_CONTROLS_NUM", length(controls.samples.all)))  
         cat("\n",controls.num, "samples are used as contorols ...\n")         
         answ <- 'y' #as.character(Sys.getenv("R_PERM_ANSW", 'y'))    
      } else if (read.from == 'stdin'){ 
         cat("\nYou have",length(cases.samples),"cases samples, which leave you",length(controls.samples.all),
             "controls. \nHow many samples you want to use for controls (the rest is use for estimating local polulation AF):")
         controls.num <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1, strip.white =T)
         cat("\nYou want to use this number of controls for all permutations (answer should be y or n):")
         answ <- scan(read.from,what=character(),n=1,quiet=T, strip.white =T)
      } else {
         if(permut ==1){
            controls.num <- length(controls.samples.all)
         }else{
            controls.num <- length(cases.samples)
         }
         answ <- 'y'
      }
      cat(paste0("Number of controls used: ", controls.num), file=logFile, append=TRUE, sep = "\n")      
      if(answ=='y'){
         same.perm <-  TRUE
      }       
   }
   
   #save.image(file='temp.stop.RData')
   annot.cols <- which(grepl('#',colnames(mutations.df)) )
   #samples.cols <- which(!grepl('#',colnames(mut.all)))
   
   
   if (controls.num < length(controls.samples.all)){     
      for (ii in pp){
         controls.samples[[ii]] <- sort(sample(controls.samples.all,controls.num))
         af.samples[[ii]] <- controls.samples.all[!controls.samples.all %in% controls.samples[[ii]]] 
         
         # - calculate local AF from mutations.df (set of all snps maybe before)
         pop.cols <- which(colnames(mutations.df) %in% af.samples[[ii]])
         pop.mut <- mutations.df[,c( annot.cols , pop.cols)]
         rownames(pop.mut) <- pop.mut$`#snpID`
         #freq <- rowSums(pop.mut[,-annot.cols], na.rm=T)/(2*length(af.samples[[ii]]))
         freq <- rowSums(pop.mut[,-annot.cols], na.rm=T)/(2*apply(pop.mut[,-annot.cols], 1 , function(x) sum(!is.na(x))  ))
         # pop.mut$`#Freq` <- apply(cbind(freq,(1-freq)),1,min) # corrected minor af
         pop.mut$`#LocalAF` <- sapply(freq, function(x) min(x, 1-x) ) # corrected minor af      
         #pop.samp.annot <- which(grepl('#',colnames(pop.mut)))
         #pop.samp.cols <- which(!grepl('#',colnames(pop.mut)))
         #pop.mut <- pop.mut[,c(pop.samp.annot, pop.samp.cols)]        
         pop.mut.df[[ii]] <- pop.mut[,c('#snpID', '#LocalAF')]
         #samples.start.pop <- sample.start+1
         rm(pop.cols)
         rm(pop.mut) 
         rm(freq)
      }     
   } else if (controls.num == length(controls.samples.all)) {
      controls.samples <- list()
      af.samples <- list()
      pop.mut.df <- list()    
      for (ii in pp){
         controls.samples[[ii]] <- controls.samples.all
         af.samples[[ii]] <- controls.samples.all     
         pop.mut.df[[ii]] <- pop.mut.df.all.contr
      }   
   } else {
      stop('You have to choose number of cases less then available number of total cases')
   }
   
   
   ########################
   #  fileter mutations 
   ########################  
   mutations.dam <- list()
   for(ii in pp){
      # Spanish frequency below af.thr
      mutations.df.fil <- merge(x=mutations.df, y=pop.mut.df[[ii]], by='#snpID', all.x=T)
      
      mut.fil.samp.annot <- which(grepl('#',colnames(mutations.df.fil)))
      mut.fil.samp.cols <- which(!grepl('#',colnames(mutations.df.fil)))
      mutations.df.fil <- mutations.df.fil[,c(mut.fil.samp.annot,mut.fil.samp.cols)]      
      mut.fil.samp.annot <- which(grepl('#',colnames(mutations.df.fil)))
      mut.fil.samp.cols <- which(!grepl('#',colnames(mutations.df.fil)))
      #mut.fil.samp.cols.cases <- which(colnames(mutations.df.fil) %in% cases.samples)
      #  mut.fil.samp.cols.contorls <- which(colnames(mutations.df.fil) %in% controls.samples[[ii]])
      #mutations.df.fil <- mutations.df.fil[,c(1:(sample.start-1), ncol(mutations.df.fil),sample.start:(ncol(mutations.df.fil)-1))]
      #sample.start.fil <- which('#LocalAF'==colnames(mutations.df.fil)) +1 
      
      mutations.df.fil <- mutations.df.fil[mutations.df.fil$`#LocalAF` < af.thr,]
      
      #exclude samples for estimating local af
      if (controls.num < length(controls.samples.all)){
         mutations.df.fil <- mutations.df.fil[,!colnames(mutations.df.fil) %in%  af.samples[[ii]] ]
      }
      #save.image('temp.2.Rdata')
      # take only rows which have at least one mutationin in cases and controls
      mut.fil.samp.cols <- which(!grepl('#',colnames(mutations.df.fil)))
      ind <- rowSums(mutations.df.fil[,mut.fil.samp.cols], na.rm=T)
      rownames(mutations.df.fil) <- mutations.df.fil$'#snpID'
      mutations.dam[[ii]] <- mutations.df.fil[ind!=0,]
      cat('In permutation', ii, 'there is', nrow(mutations.dam[[ii]]) ,' snps/indels ...\n' )
      rm(mutations.df.fil)
      #rm(mut.fil.samp.cols.cases)   
      #rm(mut.fil.samp.cols.contorls)
      rm(mut.fil.samp.annot)
      rm(mut.fil.samp.cols)
      rm(ind)
   } 
   
   
   # X matrix, same for all genes
   # SOLVE WHEN THERE IS NO COVARIATES PROVIDED ---  solved, but not well tested
   X <- list()
   for(ii in pp){
      Xi <- pca.x[c(cases.samples,controls.samples[[ii]]),, drop=F]
      Xi$intercept <- 1
      if(ncol(Xi)>2){
         Xi<- Xi[,c(ncol(Xi),2:(ncol(Xi)-1))] 
      } else {
         Xi<- Xi[,c(ncol(Xi)), drop=F]          
      }
      X[[ii]] <- Xi
      rm(Xi)
   }
   #save.image('temp.3.Rdata')
   
   ##############################################################################################
   # INLA
   ##############################################################################################      
   
   if(INLA_log){
      cat("\n Starting test inla!\n.....")
      df.perm.inla <- mclapply(pp, function(pi) perm.mc.inla(genes.all, mutations.dam[[pi]],
                                                             cases.samples, controls.samples[[pi]],
                                                             pop.mut.df[[pi]], X[[pi]], 
                                                             col.snps.numeric, col.snps.dummy, agg.col, pi), mc.cores=mc)
      
      
      
      if (p == 1 ){
         df.all.inla <- Reduce(merge, lapply(df.perm.inla, function(x) data.frame(x, genes = row.names(x))))      
         
         
         cc <- as.numeric(grep('^dic.diff_', colnames(df.all.inla) ))
         if (length(cc)>1 ){
            df.all.inla2 <- df.all.inla[order(rowMeans(df.all.inla[,cc], na.rm=T), decreasing=T),]
         } else {
            df.all.inla2 <- df.all.inla[order(df.all.inla[,cc], decreasing=T),]
         }
         
         #head(df.all.inla2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         #out.file <- "INLA"
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('INLA_logistic_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.inla2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.inla2)
         out.file.old.inla <- out.file
      } else {
         df.temp <- Reduce(merge, lapply(df.perm.inla, function(x) data.frame(x, genes = row.names(x))))  
         df.all.inla <- merge(x=df.all.inla, y=df.temp, by="genes", all=T)   
         rm(df.temp)
         
         cc <- as.numeric(grep('^cpo.diff_', colnames(df.all.inla) ))
         if (length(cc)>1 ){
            df.all.inla2 <- df.all.inla[order(rowMeans(df.all.inla[,cc], na.rm=T), decreasing=T),]
         } else {
            df.all.inla2 <- df.all.inla[order(df.all.inla[,cc], decreasing=T),]
         }
         
         #head(df.all.inla2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('INLA_logistic_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.inla2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.inla2)
         
         system(paste0('rm ',out.file.old.inla) )
         out.file.old.inla <- out.file
      }
      
      mclapply(pp, function(pi) system(paste('rm inla_perm',pi, sep='_')), mc.cores=mc)
      
      
   }
   
   
   ##############################################################################################
   # MIST
   ##############################################################################################      
   
   if(MIST_log){
      cat("\n Starting test mist!\n.....")
      df.perm.mist <- mclapply(pp, function(pi) perm.mc.mist(genes.all, mutations.dam[[pi]],
                                                             cases.samples, controls.samples[[pi]],
                                                             pop.mut.df[[pi]], X[[pi]], 
                                                             col.snps.numeric, col.snps.dummy, agg.col, pi), mc.cores=mc)
      
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
         out.file.old.mist <- out.file
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
         
         system(paste0('rm ',out.file.old.mist) )
         out.file.old.mist <- out.file
      }
      mclapply(pp, function(pi) system(paste('rm mist_perm',pi, sep='_')), mc.cores=mc)
      
   }
   
   
   ##############################################################################################
   # SKAT-O
   ##############################################################################################
   
   if(SKATO_log){
      cat("\n Starting skato!\n.....")
      df.perm.skato <- mclapply(pp, function(pi) perm.mc.skato(genes.all, mutations.dam[[pi]],
                                                               cases.samples, controls.samples[[pi]],
                                                               pop.mut.df[[pi]], X[[pi]], agg.col, pi), mc.cores=mc)
      
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
         out.file.old.skato <- out.file
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
         
         system(paste0('rm ',out.file.old.skato) )
         out.file.old.skato <- out.file
      }
      mclapply(pp, function(pi) system(paste('rm skat0_perm',pi, sep='_')), mc.cores=mc)
      
   }
   
   ##############################################################################################
   # BURDEN-O
   ##############################################################################################
   
   if(BURDEN_log){
      cat("\n Starting burden!\n.....")
      df.perm.burden <- mclapply(pp, function(pi) perm.mc.burden(genes.all, mutations.dam[[pi]],
                                                                 cases.samples, controls.samples[[pi]],
                                                                 pop.mut.df[[pi]], X[[pi]], agg.col, pi), mc.cores=mc)
      
      if (p == 1 ){
         df.all.burden <- Reduce(merge, lapply(df.perm.burden, function(x) data.frame(x, genes = row.names(x))))      
         
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.burden) ))
         if (length(cc)>1 ){
            df.all.burden2 <- df.all.burden[order(rowMeans(df.all.burden[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.burden2 <- df.all.burden[order(df.all.burden[,cc], decreasing=F),]
         }
         
         #head(df.all.burden2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         #out.file <- "BURDEN"
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('BURDEN_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.burden2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.burden2)
         out.file.old.burden <- out.file
      } else {
         df.temp <- Reduce(merge, lapply(df.perm.burden, function(x) data.frame(x, genes = row.names(x))))  
         df.all.burden <- merge(x=df.all.burden, y=df.temp, by="genes", all=T)   
         rm(df.temp)
         
         cc <- as.numeric(grep('^p.val.overall_', colnames(df.all.burden) ))
         if (length(cc)>1 ){
            df.all.burden2 <- df.all.burden[order(rowMeans(df.all.burden[,cc], na.rm=T), decreasing=F),]
         } else {
            df.all.burden2 <- df.all.burden[order(df.all.burden[,cc], decreasing=F),]
         }
         
         #head(df.all.burden2)
         out.file <- unlist(strsplit(mutation.file, '/'))
         out.file <- out.file[length(out.file)]
         
         out.file <- paste('BURDEN_results_AF',af.thr,'_Genes_',length(genes.all),'_perm',p*mc,'_',out.file,sep='')
         write.table(df.all.burden2, file=out.file, row.names=F, sep='\t', quote=F)
         rm(df.all.burden2)
         
         system(paste0('rm ',out.file.old.burden) )
         out.file.old.burden <- out.file         
      }
      mclapply(pp, function(pi) system(paste('rm burden_perm',pi, sep='_')), mc.cores=mc)
      
   }   
   
   
   ##############################################################################################
   # KBAC
   ##############################################################################################
   
   if(KBAC_log) {
      cat("\n Starting kbac!\n.....")
      df.perm.kbac <- mclapply(pp, function(pi) perm.mc.kbac(genes.all, mutations.dam[[pi]],
                                                             cases.samples, controls.samples[[pi]],
                                                             pop.mut.df[[pi]], X[[pi]], agg.col, pi), mc.cores=mc)
      
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
         
         out.file.old.kbac <- out.file
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
         
         system(paste0('rm ',out.file.old.kbac) )
         out.file.old.kbac <- out.file
      }
      mclapply(pp, function(pi) system(paste('rm kbac_perm',pi, sep='_')), mc.cores=mc)
      
   }
   
   
}

save(af.samples, file='af.samples.Rdata')
save(controls.samples, file='controls.samples.Rdata')

cat(paste0("In file af.samples.Rdata are saved (R object) af samples IDs used for AF estimation"), file=logFile, append=TRUE, sep = "\n")
cat(paste0("In file controls.samples.Rdata are saved (R object) af samples IDs used as controls"), file=logFile, append=TRUE, sep = "\n")
cat(paste0("Everything done!"), file=logFile, append=TRUE, sep = "\n")
