######################################################################################
# CRG 
# Authors: Hana SUSAK, Georgia ESCARAMIS 
# date: 25/06/2018
#------------------------------------------------------------------------------------
# QC check for multisample call (PCA, filtering on sample info file, #mut per sample)
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Calling the script: Rscript QC_and_filtering_multisaple_call_v2.R --h
#------------------------------------------------------------------------------------
######################################################################################

rm(list=ls())
time.f <- Sys.time()
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(argparse))
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
#suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser()

parser$add_argument("-s", "--snp_file", type="character", help="input snp file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-i", "--indel_file", type="character", help="input indel file", metavar="file", nargs=1)
parser$add_argument("-d", "--sample_desc_file", type="character", help="samples description file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-c", "--correct_sample_names", type="character", help="file with corrected names", metavar="file", nargs=1)

read.from <- 'stdin'
# read.from <- ''

# setting log file and output folder
out_foleder <-  paste0('QC_results_folder_',format(time.f, format = "%Y-%m-%d-%H%M%S") )

#gc()
## get command line arguments
args <- parser$parse_args()
snp.file.name <- args$snp_file
indel.file.name <- args$indel_file
samples.info.file <- args$sample_desc_file
sample.correct.file.name <- args$correct_sample_names
# snp.file.name <- '/no_backup/GD/projects/HumanDisease/HaplotypeCaller_march15/genotype/genotype_220515/NEW_FEBRUARY_FS/intersection/agi50_agi71_nimv3/snps/old_noXY/snps_agi50_agi71_nimv3_parsed.csv'
# indel.file.name <- '/no_backup/GD/projects/HumanDisease/HaplotypeCaller_march15/genotype/genotype_220515/NEW_FEBRUARY_FS/intersection/agi50_agi71_nimv3/indels/old_noXY/indels_agi50_agi71_nimv3_parsed.csv'
# samples.info.file <- '/no_backup/GD/projects/CLL/germ_analysis/data/Samples_info_all_clean.5_fin.txt'
# sample.correct.file.name <- NULL

# snp.file.name <- '/no_backup/GD/projects/RVAS/1kg_analysis/data/1kg.all.clean.merged.txt'
# samples.info.file <-'/no_backup/GD/projects/RVAS/1kg_analysis/data/1kg_info_parsed_rel.txt'

# snp.file.name <- '/no_backup/GD/projects/RVAS/data/toy_example_snps_1000_samples.txt'
# indel.file.name <- '/no_backup/GD/projects/RVAS/data/toy_example_indels_1000_samples.txt'
# samples.info.file <- '/no_backup/GD/projects/RVAS/data/Samples_info_all_clean.5_fin.txt'
# sample.correct.file.name <- NULL

# snp.file.name <- '/users/xe/sozkan/parsing/snps/snps_agi50_agi71_nimv3_parsed_headerok.csv'
# indel.file.name <- '/users/xe/sozkan/parsing/indels/indels_agi50_agi71_nimv3_parsed_headerok.csv'
# samples.info.file <- '/no_backup/GD/projects/CLL/germ_analysis/data/Samples_info_all_clean.5_fin.txt'
# sample.correct.file.name <- NULL

######################################################################################
# read input
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n ----------------------------------- Reading input files -------------------------------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# normalize paths 
snp.file.name <- normalizePath(snp.file.name)
samples.info.file <- normalizePath(samples.info.file)
if(!is.null(indel.file.name)){
   indel.file.name <- normalizePath(indel.file.name)   
}
if(!is.null(sample.correct.file.name)){
   sample.correct.file.name <- normalizePath(sample.correct.file.name)   
}

dir.create(out_foleder)
setwd(out_foleder)
logFile <-  paste0('QC_analysis_log_file_',format(time.f, format = "%Y-%m-%d-%H%M%S") , '.txt')
cat(timestamp(quiet=T), file=logFile, append=FALSE, sep = "\n")
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('------------------------------------ Reading input files -------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

cat("\nAll results you can find in folder:", getwd(),"\n" )

cat(paste('Input SNPs file:', snp.file.name, sep='\t'), file=logFile, append=TRUE, sep = "\n")
cat(paste('Input InDels file:', indel.file.name, sep='\t'), file=logFile, append=TRUE, sep = "\n")
cat(paste('Input Samples info file:', samples.info.file, sep='\t'), file=logFile, append=TRUE, sep = "\n")
cat(paste('Input correct sample names file:', sample.correct.file.name, sep='\t'), file=logFile, append=TRUE, sep = "\n")

if(!is.null(sample.correct.file.name)){
   correct.names <- read.table(file=sample.correct.file.name, header=T, comment.char="", quote="", sep='\t', stringsAsFactors=F, check.names=F)
   row.names(correct.names) <- correct.names$wrong 
}

#till here 
cat("start! \n .....")


######################################################################################
# read snps and indels
######################################################################################
# SNPs
mut.all <- read.table(file=snp.file.name, header=T, comment.char="", quote="", sep='\t',  
                       # nrows=2000000,
                       stringsAsFactors=F, check.names=F)
cat("\nIn the input file there are", nrow(mut.all), "unique snps  \n ....." )

refseq.pos <- which(grepl('\\(Refseq\\)',colnames(mut.all)))
new.names <- gsub('\\(Refseq\\)','',colnames(mut.all)[refseq.pos])
colnames(mut.all)[refseq.pos] <- new.names

rm(new.names)
rm(refseq.pos)
################################################
# select columns user wants + obligatory 4
################################################
annot.cols <-  which(grepl('#',colnames(mut.all)))
#samples.cols <- which(!grepl('#',colnames(mut.all)))

cat("\n\nWhich columns you want to keep (later to fileter on or aggregate) :", 
    paste('\n',1:length(annot.cols),'.',colnames(mut.all)[annot.cols], sep=''),
    "\nWrite numbers that are in front of the columns you want to keep (space separated) : ")
col.keep <- scan(read.from,what=numeric(),quiet=TRUE, nlines=1)
col.keep <- unique(col.keep)
col.keep <- col.keep[col.keep %in% (annot.cols)]

annot.cols.names <- colnames(mut.all)[annot.cols]
annot.cols.names <- annot.cols.names[col.keep]
   
annot.names.keep.must <- c('#Chr','#Position','#Reference','#Alteration')
annot.names.keep <- unique(c(annot.names.keep.must,annot.cols.names))

if(any(!annot.names.keep %in% colnames(mut.all)) ){
   message(paste('You should have this column names in snps/indels file:', paste(annot.names.keep, collapse=', ')))
   stop('There are missing or wrongly named columns in your file')
} else {
   message('Annotation columns to be kept are: ',paste(annot.names.keep, collapse=', '))
}

rm(annot.cols.names)
rm(col.keep)
rm(annot.names.keep.must)

#log file of kept columns
cat(paste('Columns you kept for downstream analyis: ', paste0(annot.names.keep, collapse=", "), sep=' '), file=logFile, append=TRUE, sep = "\n")

annot.cols <- match(annot.names.keep, colnames(mut.all))
samples.cols <- which(!grepl('#',colnames(mut.all)))

mut.all <- mut.all[,c(annot.cols, samples.cols)]

if(!is.null(indel.file.name)){
   # InDels
   indels.all <- read.table(file=indel.file.name, header=T, comment.char="", quote="", sep='\t',  # nrows=100000, 
                            stringsAsFactors=F, check.names=F)
   cat("\nIn input file there are", nrow(indels.all), "unique indels \n ....." )
   
   refseq.pos <- which(grepl('\\(Refseq\\)',colnames(indels.all)))
   new.names <- gsub('\\(Refseq\\)','',colnames(indels.all)[refseq.pos])
   colnames(indels.all)[refseq.pos] <- new.names
   
   rm(refseq.pos)
   rm(new.names)
   
   if(any(!annot.names.keep %in% colnames(indels.all)) ){
      message(paste('You should have this column names in snps/indels file:', paste(annot.names.keep, collapse=', ')))
      stop('There are missing or wrongly named columns in your file')
   }
   
   annot.cols <- match(annot.names.keep, colnames(indels.all))
   samples.cols <- which(!grepl('#',colnames(indels.all)))
   
   indels.all <- indels.all[,c(annot.cols, samples.cols)]
   
   #merge data sets
   indels.all$`#type` <- 'indel'
   mut.all$`#type` <- 'snp'
   
   if (ncol(indels.all)==ncol(mut.all) & all(colnames(mut.all) == colnames(indels.all))){
      mut.all <- rbind(mut.all, indels.all)
      
      annot.cols <- which(grepl('#',colnames(mut.all)))
      samples.cols <- which(!grepl('#',colnames(mut.all)))     
      mut.all <- mut.all[,c(annot.cols, samples.cols)]
      
      annot.cols <- which(grepl('#',colnames(mut.all)))
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      
      #rm(snps.all)
      rm(indels.all) 
   } else {
      stop('Column names in snps and indels file are not matching!')
   } 
} else {
      mut.all$`#type` <- 'snp'
   #mut.all <- snps.all
   annot.cols <- which(grepl('#',colnames(mut.all)))
   samples.cols <- which(!grepl('#',colnames(mut.all)))
   
   cat("\n\nDo you have column in snp file which representd annotation for snps vs indels (answer y or n)?")
   snpVSindel.col <- scan(read.from,what=character(),n=1,quiet=TRUE)
   if(snpVSindel.col=='y'){
      cols.types <- sapply(mut.all[,annot.cols], class)
      char.cols <- cols.types[cols.types =='character']
      cat("\n\nWhich columns have identifyer for snps or idenls  :", 
          paste('\n',1:length(char.cols),'.',names(char.cols), sep=''),
          "\nWrite number that is in front of the column you want to choose : ")
      col.indel.snps <- scan(read.from,what=numeric(),n=1, quiet=TRUE, nlines=1)
      indel.snps <- names(char.cols)[col.indel.snps]
      options <- as.character(unique(mut.all[,indel.snps]))
      options[is.na(options)] <- 'NA_Missing_info'
      cat("\n\nIn column",indel.snps,"there are folowing fields:")
      counts.tb <- table(mut.all[,indel.snps], useNA='always')
      names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
      counts.tb <- counts.tb[options]
      cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' samples)','\n', sep=''))
      cat("\nEnter the numbers which are in front of the options you consider as indels (space separated): ")
      indel.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
      indel.anot.names <- options[indel.options]
  
      if( 'NA_Missing_info' %in% options[indel.options]){
         mut.all[ is.na(mut.all[,indel.snps]), '#type'] <-  'indel'
      }
      mut.all[mut.all[,indel.snps] %in% indel.anot.names &
                 !is.na(mut.all[,indel.snps]), '#type'] <-  'indel'
      
      #log file - indels 
      cat(paste('Mutations marked as',paste0(indel.anot.names, collapse=', '), 'in column',indel.snps, 'are maked as indels' , sep=' '), file=logFile, append=TRUE, sep = "\n")
      
      rm(options)
      rm(col.indel.snps)
      rm(indel.snps)
      rm(counts.tb)
      rm(indel.options)
      rm(snpVSindel.col)
      rm(indel.anot.names)
      rm(char.cols)
      rm(cols.types)
      
      annot.cols <- which(grepl('#',colnames(mut.all)))
      samples.cols <- which(!grepl('#',colnames(mut.all)))   
     
   } 
   mut.all <- mut.all[,c(annot.cols, samples.cols)]
   annot.cols <- which(grepl('#',colnames(mut.all)))
   samples.cols <- which(!grepl('#',colnames(mut.all)))
   #rm(snps.all)   
}
rm(annot.names.keep)

######################################################################################
# filtering SNPs/InDels based on annotation
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n --------------------------------- Filtering by SNPs/InDels ----------------------------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('---------------------------------- Filtering by SNPs/InDels ----------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

cat("\nDo you want to filter by any columns of the mutation file [ by muation annotation, not by samples ](answer y or n)? \nAnswer: ")
filt.ind <- scan(read.from,what=character(),n=1,quiet=TRUE, nlines=1)
filt.ind <- tolower(substr(filt.ind,1,1))
if(filt.ind =='y'){
   annot.cols <- which(grepl('#',colnames(mut.all)))
   samples.cols <- which(!grepl('#',colnames(mut.all)))
   
   cols.types <- sapply(mut.all[,annot.cols], class)
   char.cols <- cols.types[cols.types =='character']
   num.cols <- cols.types[cols.types =='numeric' |cols.types =='integer' ]
   rm(cols.types)
   
   #character columns
   cat("\nYou can choose chaarachter  colums to filter by:", 
       paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
   cat("\nEnter the column numbers, space separated: ")
   x <- scan(read.from,what=integer(),nmax=length(char.cols),quiet=TRUE, nlines=1)
   x <- unique(x)
   x <- x[x<= length(char.cols)]       
   # log file separator
   cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
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
         options <- as.character(unique(mut.all[,temp.col]))
         options[is.na(options)] <- 'NA_Missing_info'
         cat("\n\nIn column",names(char.cols)[i],"there are folowing fields:")
         counts.tb <- table(mut.all[,temp.col], useNA='always')
         names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
         counts.tb <- counts.tb[options]
         cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' mutations)','\n', sep=''))
         cat("Enter the numbers which are in front of the options you DON'T want to keep (space separated): ")
         exclude.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
         exclude.options <- unique(exclude.options)
         exclude.options <- exclude.options[exclude.options <= length(options)]
         if (length(exclude.options) > 0) {
            if( 'NA_Missing_info' %in% options[exclude.options]){
               mut.all <-  mut.all[!is.na(as.character(mut.all[,temp.col])) , ]
            }
            mut.all <- mut.all[! mut.all[,temp.col] %in% options[exclude.options] | is.na( mut.all[,temp.col]),  ]
            cat("\nMutations whos", temp.col,"is in:", paste(options[exclude.options], collapse=', '),"- will be excluded")  
            #log file - snps character annot keept
            cat(paste("Mutations are keept if",temp.col,"is matching:",paste(options[-exclude.options], collapse=', ' ), sep=' '), 
                 file=logFile, append=TRUE, sep = "\n")            
         } else {
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
   
   # log file separator
   cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")

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
         quant.temp.col <- quantile(mut.all[,temp.col], na.rm=T)
         NA.temp.col <- sum(is.na(mut.all[,temp.col]))
         if(NA.temp.col >0 ){
           cat("\nIn column",temp.col,"therea are",NA.temp.col,'Na values')
           cat("\n\nDo you want to filter out these mutations (answer y or n)? \nAnswer:")
           NA.ind <- scan(read.from,what=character(),n=1,quiet=TRUE)
           if (NA.ind=='y'){
              mut.all <- mut.all[!is.na(mut.all[,temp.col]),  ]
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
             mut.all <- mut.all[mut.all[,temp.col] >= trehsold | is.na(mut.all[,temp.col]),  ]
             message('Everything lower then ',trehsold,' in column ', temp.col, ' will be removed')
             #log file - snps numeric annot keept
             cat(paste("Mutations are keept if",temp.col,"is greater or equal then:", trehsold, sep=' '), 
                 file=logFile, append=TRUE, sep = "\n")    
         } else if (num.ind == 'h') {
             mut.all <- mut.all[mut.all[,temp.col] <= trehsold | is.na(mut.all[,temp.col]),  ]
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
}       
rm(filt.ind)

annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))
# remove samples with no called mutations
remove.sam <- which(colSums(mut.all[,samples.cols],na.rm=T) == 0)
if(length(remove.sam)>0){
   cat('There are samples with 0 mutations:',paste(names(remove.sam), collapse=', '),'\nDo you want to remove them (y or n): ')
   answ.zero <- scan(read.from,what=character(),n=1,quiet=TRUE)
   answ.zero <- tolower(substr(answ.zero,1,1))
   if(answ.zero=='y'){
      cols.to.exclude <- which(colnames(mut.all) %in% names(remove.sam)) 
      mut.all <- mut.all[,-cols.to.exclude] 
      cat('They are now removed.\n')
      #log file - snps numeric annot keept
      cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
      cat(paste("Samples with 0 mutations are excluded. These samples are:", paste0(names(remove.sam), collapse=', '), sep=' '), 
          file=logFile, append=TRUE, sep = "\n") 
   }
   annot.cols <- which(grepl('#',colnames(mut.all)) )
   samples.cols <- which(!grepl('#',colnames(mut.all)))
   rm(answ.zero)
}
rm(remove.sam)

######################################################################################
# filtering snps based on NA precentege
######################################################################################

na.prec.per.snp <- apply(mut.all[,samples.cols], 1, function(x) sum(is.na(x)/length(x)))
quant.snps.na <- quantile(na.prec.per.snp, probs=c(0,0.25, 0.5, 0.75, 0.9, 0.95, 1))
cat("\nNAs [fractions] per mutations quantiles are :", 
    paste(paste0('\n',quant.snps.na), names(quant.snps.na), sep=' at qunatile ' ) )
cat("\n\nDo you want to filter out mutations with some fractoin of NA values (y or n)? \nAnswer: ")
na.snp.ind <- scan(read.from,what=character(),n=1,quiet=TRUE)
na.snp.ind <- tolower(na.snp.ind)
na.snp.ind <- substr(na.snp.ind,1,1)
if(na.snp.ind=='y'){
   cat("\nWhat is treshold you want to use (in fraction). Everythig about   : ")
   trehsold <- scan(read.from,what=numeric(),n=1,quiet=TRUE)      
   tot.snps <- nrow(mut.all)
   rm.na.cols <- which(na.prec.per.snp >= trehsold)
   mut.all <- mut.all[-rm.na.cols,]
   message('Every mutation which have more then ',trehsold*100,'% of NA calls in samples will be filtered out')
   message('There are  ',length(rm.na.cols),' number ( ', round((length(rm.na.cols)/tot.snps)*100, digits=2), '% ) of such mutations.') #log file - snps numeric annot keept   
   cat(paste("Every mutation which have more then",trehsold*100, "% of NA calls in samples will be filtered out.", sep=' '),file=logFile, append=TRUE, sep = " ")     
   cat(paste('There are ',length(rm.na.cols),'(', round((length(rm.na.cols)/tot.snps)*100, digits=2),
                 '% ) of such mutations.', sep=' '), file=logFile, append=TRUE, sep = "\n")   
   rm(trehsold)
   rm(rm.na.cols)
   rm(tot.snps)
}
 
rm(na.snp.ind)
rm(quant.snps.na)
rm(na.prec.per.snp)
######################################################################################
# correcting some names
######################################################################################
if(exists('correct.names')){
   cat("\n ---------------------------------------------------------------------------------------------\n")
   cat("\n------------------------------------- Correcting names  --------------------------------------\n")
   
   # lof file
   cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
   cat('-------------------------------------- Correcting names --------------------------------------', file=logFile, append=TRUE, sep = "\n")
   cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")
   
   # correct some snps names
   snps.all.colnames <- colnames(mut.all)
   cor.index <- which(snps.all.colnames %in% rownames(correct.names))
   cor.name <- correct.names[ snps.all.colnames[snps.all.colnames %in% rownames(correct.names)] , 2]
   
   if (length(cor.index) > 0){
      cat("\nSample names in files: ", paste(colnames(mut.all)[cor.index], collapse=', '), 
          '\nare corrected to:',paste(cor.name, collapse=', ')," \n .....")
      colnames(mut.all)[cor.index] <- cor.name  
      # log file
      cat(paste("Sample names in files:", paste(colnames(mut.all)[cor.index], collapse=', '),
                "\nare corrected to:",paste(cor.name, collapse=', '), sep=' '), file=logFile, append=TRUE, sep = "\n") 
   }
   rm(cor.index)
   rm(cor.name)
   rm(snps.all.colnames)
   rm(correct.names)
   cat("\n------------------------------ Done with correcting names ------------------------------------\n")
   cat("\n ---------------------------------------------------------------------------------------------\n")
} else {
   rm(sample.correct.file.name)
}


#save.image('temp1.RData')
#load('temp1.RData')
######################################################################################
# filtering by Samples info file - removing samples
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n---------------------------------- Filtering by SAMPLE IDs -----------------------------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('---------------------------------- Filtering by SAMPLE IDs -----------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

samples.to.exclude <- c()

# collect  info for samples
if(exists('samples.info.file')){
   samples.info <- read.table(file=samples.info.file, header=T, sep="\t", quote="", fill=T)
   message("Samples info/description file must contain samples IDs")
   cat("What is column with samples id (like in multisample snps/indels file) :", 
       paste('\n',1:ncol(samples.info),'.',colnames(samples.info), sep=''),
       "\nChoose number: ")
   answ.col <- scan(read.from,what=integer(),n=1,quiet=TRUE)
   answ.col <- answ.col[answ.col<= ncol(samples.info)]       
   answ.name <- colnames(samples.info)[answ.col]
   cat("\nYou chose that samples id column is:", answ.name, ". ") 
   cat(paste("You chose that samples id column is:", answ.name),  file=logFile, append=TRUE, sep = "\n")
   cat("\n ---------------------------------------------------------------------------------------------\n")
   answ <- 'y'
   #colnames(samples.info)[answ.col] <- 'samples'
   if( any(!colnames(mut.all)[samples.cols] %in% samples.info[,answ.col])){
      message(' In samples info file there are less samples then  samples from mutation file.\n Mutation file will be subselected to interesect with the samples from info file')   
      annot.cols <- which(grepl('#',colnames(mut.all)) )
      keep.samp.info <- which(colnames(mut.all) %in% samples.info[,answ.name])      
      mut.all <- mut.all[,c(annot.cols,keep.samp.info )]
      # recalculate columns again
      annot.cols <- which(grepl('#',colnames(mut.all)) )
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      rm(keep.samp.info)
      
      # remove snps with no mutation
      remove.snps <- which(rowSums(mut.all[,samples.cols], na.rm=T) == 0)
      if(length(remove.snps)>0){
         mut.all <- mut.all[-remove.snps,] 
      }
      rm(remove.snps) 
      
   } 
   samples.info <- samples.info[samples.info[,answ.name] %in%  colnames(mut.all)[samples.cols],]      
   cat("\n\nYou are starting with ", nrow(samples.info), " samples")
   cat(paste("You are starting with:", nrow(samples.info), " samples"),  file=logFile, append=TRUE, sep = "\n")
   
   
   cat("\n\nYou can filter by this colums :", 
       paste('\n',1:ncol(samples.info),'.',colnames(samples.info), sep=''),
       "\nYou want to filer (y or n): ")
   answ <- scan(read.from,what=character(),n=1,quiet=TRUE)
   answ <- tolower(substr(answ,1,1))   
   if(answ=='y'){
      cat("\nEnter the column numbers, space separated: ")
      x <- scan(read.from,what=integer(),nmax=ncol(samples.info),quiet=TRUE, nlines=1)
      x <- unique(x)
      x <- x[x<= ncol(samples.info)]       
      if (length(x) <= 0){
         cat("\nYou or chose not to filter or wrong columns \n.....")
      } else {
         cat("\nColumn(s) you choose are:",paste(colnames(samples.info)[x], collapse=', ' ), "\n.....")
         # log file
         cat(paste("Column(s) you choose are:", paste(colnames(samples.info)[x], collapse=', ' )),  file=logFile, append=TRUE, sep = "\n")
         for (i in x){
            options <- as.character(unique(samples.info[,i]))
            options[is.na(options)] <- 'NA_Missing_info'
            cat("\n\nIn column",colnames(samples.info)[i],"there are folowing fields:")
            counts.tb <- table(samples.info[,i], useNA='always')
            names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
            counts.tb <- counts.tb[options]
            cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' samples)','\n', sep=''))
            cat("\nEnter the numbers which are in front of the options you DON'T want to keep (space separated): ")
            exclude.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
            exclude.options <- unique(exclude.options)
            exclude.options <- exclude.options[exclude.options <= length(options)]
            if (length(exclude.options) > 0) {
               if( 'NA_Missing_info' %in% options[exclude.options]){
                  samples.to.exclude.temp <-  as.character(samples.info[is.na(as.character(samples.info[,i])) ,answ.name ])
                  samples.to.exclude <- unique(c(samples.to.exclude, samples.to.exclude.temp)) 
               }
               samples.to.exclude.temp <-  as.character(samples.info[as.character(samples.info[,i]) %in% options[exclude.options],answ.name ])
               samples.to.exclude <- unique(c(samples.to.exclude, samples.to.exclude.temp)) 
               cat("\nSamples whos", colnames(samples.info)[i],
                        "is in:", paste(options[exclude.options], collapse=', '),
                        "- will be excluded") 
               # log file
               cat(paste("Samples that are kept have", colnames(samples.info)[i], "matching:" ,
                         paste(options[-exclude.options], collapse=', ')),  file=logFile, append=TRUE, sep = "\n")                  
            }  else {
               cat("\n Nothing from column",colnames(samples.info)[i],"will be excluded ")
               cat(paste("Nothing from column",colnames(samples.info)[i],"will be excluded"),  file=logFile, append=TRUE, sep = "\n") 
            }     
            cat("\n ---------------------------------------------------------------------------------------------")
         }
         rm(samples.to.exclude.temp)
         rm(exclude.options)         
         rm(options)
         rm(counts.tb)
         rm(i)
      }  
      rm(x)
   } 
   rm(answ)
   rm(answ.col)
} else {
   stop('You have to provide Samples info file')
}

cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
cat("\nTotal number of samples to exclude is",length(cols.to.exclude)," \n.....")

if (length(cols.to.exclude) > 0){
   mut.all <- mut.all[,-cols.to.exclude]
   samples.info2 <- samples.info[!samples.info[,answ.name] %in% samples.to.exclude,]
   samples.info2 <- droplevels(samples.info2)    
} else {
   samples.info2 <- samples.info
}
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

remove.snps <- which(rowSums(mut.all[,samples.cols], na.rm=T) == 0)
if(length(remove.snps)>0){
   mut.all <- mut.all[-remove.snps,] 
   cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
   cat(paste("Mutations with 0 samples mutated (or only NA genotypes) are excluded. Amount of such mutations is:",length(remove.snps), sep=' '), 
       file=logFile, append=TRUE, sep = "\n") 
}
rm(remove.snps) 
rm(samples.info)
rm(snp.file.name)
rm(indel.file.name)
rm(samples.info.file)
rm(cols.to.exclude)
#log file
cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of samples left after filtering:",nrow(samples.info2)),  file=logFile, append=TRUE, sep = "\n")
#save.image('Temp.data2.RData')

######################################################################################
# Filter outliers by number of mutations
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n --------------------- Filtering outlier samples by number of mutations ----------------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('--------------------- Filtering outlier samples by number of mutations -----------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

mut.all$`#snpID` <- apply(mut.all[,c("#Chr","#Position","#Reference","#Alteration")], 1, function(x) paste(trimws(x),collapse=';'))
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

# reorder, first annotation, then samples columns
mut.all <- mut.all[,c(annot.cols, samples.cols )]
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

cols.types <- sapply(mut.all[,annot.cols], class)
char.cols <- cols.types[cols.types =='character']
rm(cols.types)
cat("\n\nYou can choose now column to be used later for coloring mutation types in barplot:", 
    paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
cat("\nEnter the column number (or if you don't want to do this, press enter): ")
x.col <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
x.col <- unique(x.col)
x.col <- x.col[x.col<= length(char.cols)]
if (length(x.col) <= 0){
   cat("\nYou choose not to color by any annotation column \n.....")
   cc <- which(colnames(mut.all)%in% c( "#snpID"))
} else {
   col.snps.nar <- names(char.cols[x.col])
   cat("\nColumn you choose is:",paste(col.snps.nar, collapse=', ' ), "\n.....")
   cat(paste('\nColumn you choose to color by in barplot is:',paste(col.snps.nar, collapse=', ' )), file=logFile, append=TRUE, sep = "\n")
   cc <- which(colnames(mut.all)%in% c(c( "#snpID"), col.snps.nar))
}
rm(char.cols)
rm(x.col)

for(i in 1:ceiling(nrow(mut.all)/100000)){
   max.row <- min((i*100000),nrow(mut.all)) 
   #print(max.row)
   if (exists('col.snps.nar')){
      melt.mut.all.temp <- melt(mut.all[((i-1)*100000+1):max.row,c(cc,samples.cols)], id.vars=c("#snpID", col.snps.nar))
      melt.mut.all.temp <- melt.mut.all.temp[melt.mut.all.temp$value != 0 & !is.na(melt.mut.all.temp$value),]
   } else {
      melt.mut.all.temp <- melt(mut.all[((i-1)*100000+1):max.row,c(cc,samples.cols)], id.vars=c("#snpID"))
      melt.mut.all.temp <- melt.mut.all.temp[melt.mut.all.temp$value != 0 & !is.na(melt.mut.all.temp$value),]   
   }
   if (i ==1){
      melt.mut.all <- melt.mut.all.temp
      rm(melt.mut.all.temp)
      #print('added!')
   } else {
      melt.mut.all <- rbind(melt.mut.all, melt.mut.all.temp)
      rm(melt.mut.all.temp)
      #print('added!')
   }
}
rm(max.row)
rm(i)
rm(cc)
#melt.mut.all <- melt.mut.all[!is.na(melt.mut.all$value),]

all.samples <- as.character(colnames(mut.all[,samples.cols]))
melt.mut.all$variable <- factor(melt.mut.all$variable, levels=all.samples)
#save.image('Temp.data3.RData')

#########
# here plot histogram
total.per.sample <- aggregate(`#snpID` ~ variable,data=melt.mut.all, length) 
colnames(total.per.sample) <- c('samples', 'Num.Mutations')
total.per.sample$samples <- factor(total.per.sample$samples, levels=all.samples)
rm(all.samples)

total.per.sample <- total.per.sample[order(total.per.sample$Num.Mutations, decreasing=T), ]
zero.samples <- levels(total.per.sample$samples)[!levels(total.per.sample$samples) %in% total.per.sample$samples] 
if(length(zero.samples)>0){
   tot.temp <- data.frame(samples=zero.samples, Num.Mutations=0 )
   total.per.sample <- rbind(total.per.sample, tot.temp)
   rm(tot.temp)
}
rm(zero.samples)


samples.info2 <- merge(x=samples.info2, y=total.per.sample, by.x=answ.name, by.y='samples', all.x=TRUE)
rm(total.per.sample)

X11(width=14, height=8)
p1 <- ggplot(samples.info2, aes(x=Num.Mutations)) + geom_histogram() + theme_bw() + ggtitle('Number of mutations per sample histogram')
p1

cat("\nEnter minimum number of mutations per sample: ")
min.mut <- scan(read.from,what=numeric(),nmax=1,quiet=TRUE, nlines=1)
cat("\nEnter maximum number of mutations per sample: ")
max.mut <- scan(read.from,what=numeric(),nmax=1,quiet=TRUE, nlines=1)
#log file
cat(paste("Minimum number of mutations per sample:",min.mut),  file=logFile, append=TRUE, sep = "\n")
cat(paste("Maximum number of mutations per sample:",max.mut),  file=logFile, append=TRUE, sep = "\n")

dev.off()
#plot hist of per saple number of mutations
pdf(file='hist.number.mutations_old.pdf', width=15,height=12,useDingbats=F)
p1
dev.off()
rm(p1)
#log file
cat(paste("hist.number.mutations_old.pdf done! - histogram for number of mutations per sample - not final"),
    file=logFile, append=TRUE, sep = "\n")

# find samples which don't fit the criteria
samples.to.exclude.temp <- c()
if (length(min.mut) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(samples.info2[samples.info2$Num.Mutations < min.mut,answ.name]))
}
if (length(max.mut) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(samples.info2[samples.info2$Num.Mutations > max.mut,answ.name]))
}
samples.to.exclude.temp <- unique(samples.to.exclude.temp)

# remove such samples if they exist
if(length(samples.to.exclude.temp)>0){
   cat("Samples you decided to exclude because of number of mutations are:",
       paste(samples.to.exclude.temp, collapse=", "),"\n.....")
   #log file
   cat(paste("Samples exclude because of number of mutations are:", paste(samples.to.exclude.temp, collapse=", ")),
       file=logFile, append=TRUE, sep = "\n") 
   samples.to.exclude <- c(samples.to.exclude, samples.to.exclude.temp)
   rm(max.mut)
   rm(min.mut)
   
   cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
   cat("\nTotal number of samples to exclude because of number of mutations is",length(samples.to.exclude.temp),"\n.....")
   #log file
   cat(paste("Total number of samples to exclude because of number of mutations is: ",length(samples.to.exclude.temp)),
       file=logFile, append=TRUE, sep = "\n") 
   rm(samples.to.exclude.temp)
   
   if (length(cols.to.exclude) > 0){
      mut.all <- mut.all[,-cols.to.exclude]   
      melt.mut.all <- melt.mut.all[!melt.mut.all$variable %in% samples.to.exclude,]
      melt.mut.all <- droplevels(melt.mut.all)
      samples.info2 <- samples.info2[!samples.info2[,answ.name] %in% samples.to.exclude,]
      samples.info2 <- droplevels(samples.info2)
      #all.samples <- all.samples[!all.samples %in% samples.to.exclude]
      
      #recaluclate annotation and sample columns
      annot.cols <- which(grepl('#',colnames(mut.all)) )
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      
      remove.snps <- which(rowSums(mut.all[,samples.cols], na.rm=T) == 0)
      if(length(remove.snps)>0){
         mut.all <- mut.all[-remove.snps,]  
         cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
         cat(paste("Mutations with 0 (or only NAs) samples mutated are excluded. Amount of such mutations is:",length(remove.snps), sep=' '), 
             file=logFile, append=TRUE, sep = "\n") 
      }
      rm(remove.snps) 
   }
   rm(cols.to.exclude) 
} else {
   rm(samples.to.exclude.temp)
   rm(max.mut)
   rm(min.mut)
   cat("You choose not to remove any sample. \n.....")
}

#log file
cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of samples left after filtering:",nrow(samples.info2)),  file=logFile, append=TRUE, sep = "\n")

######################################################################################
# Tv/Ts ratio by sample
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n ------------------------- Filtering outlier samples by TR/TV ratio --------------------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('-------------------------- Filtering outlier samples by TR/TV ratio --------------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

melt.mut.all$tvtr <- NA
melt.mut.all[grepl(';A;G$',melt.mut.all$`#snpID`),'tvtr'] <- 'tr'
melt.mut.all[grepl(';G;A$',melt.mut.all$`#snpID`),'tvtr'] <- 'tr'
melt.mut.all[grepl(';C;T$',melt.mut.all$`#snpID`),'tvtr'] <- 'tr'
melt.mut.all[grepl(';T;C$',melt.mut.all$`#snpID`),'tvtr'] <- 'tr'

melt.mut.all[grepl(';T;G$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';G;T$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';C;A$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';A;C$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';C;G$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';G;C$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';T;A$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'
melt.mut.all[grepl(';A;T$',melt.mut.all$`#snpID`),'tvtr'] <- 'tv'

tvtr.ratio <- aggregate(`#snpID` ~ variable + tvtr,data=melt.mut.all, length) 
tvtr.ratio <- tvtr.ratio[!is.na(tvtr.ratio$tvtr) ,]
colnames(tvtr.ratio) <- c('Sample','tvtr','count')
tvtr.ratio <- dcast(tvtr.ratio,Sample ~ tvtr, fun.aggregate=sum, value.var='count' )
tvtr.ratio$tvtr_ratio <- tvtr.ratio$tr/tvtr.ratio$tv
tvtr.ratio <- tvtr.ratio[,c('Sample', 'tvtr_ratio')]
samples.info2 <- merge(x=samples.info2, y=tvtr.ratio, by.x=answ.name, by.y='Sample', all.x=TRUE)
rm(tvtr.ratio)

X11(width=14, height=8)
p1 <- ggplot(samples.info2, aes(x=tvtr_ratio)) + geom_histogram() + theme_bw() + ggtitle('Tr/Tv ratio per sample histogram')
p1

cat("\nEnter minimum tv/tr ratio (per sample): ")
min.tvtr <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter maximum tv/tr ratio (per sample): ")
max.tvtr <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
#log file
cat(paste("Minimum tv/tr ratio per sample:",min.tvtr),  file=logFile, append=TRUE, sep = "\n")
cat(paste("Maximum tv/tr ratio per sample:",max.tvtr),  file=logFile, append=TRUE, sep = "\n")

dev.off()
#plot hist of per saple number of mutations
pdf(file='hist.tvtr.ratio_old.pdf', width=15,height=12,useDingbats=F)
p1
dev.off()
rm(p1)
cat(paste("hist.tvtr.ratio_old.pdf done! - histogram for TR-TV ratio per sample - not final"),
    file=logFile, append=TRUE, sep = "\n")

# find samples with don't fit the criteria
samples.to.exclude.temp <- c()
if (length(min.tvtr) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(samples.info2[samples.info2$tvtr_ratio < min.tvtr,answ.name]))
}
if (length(max.tvtr) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(samples.info2[samples.info2$tvtr_ratio > max.tvtr,answ.name]))
}
samples.to.exclude.temp <- unique(samples.to.exclude.temp)

# remove such samples if they exist
if(length(samples.to.exclude.temp)>0){
   cat("Samples you decided to exclude because of TV-TR ratio are:",
       paste(samples.to.exclude.temp, collapse=", "),"\n.....")
   #log file
   cat(paste("Samples exclude because of TV-TR ratio are:", paste(samples.to.exclude.temp, collapse=", ")),
       file=logFile, append=TRUE, sep = "\n") 
   samples.to.exclude <- c(samples.to.exclude, samples.to.exclude.temp)
   rm(samples.to.exclude.temp)
   rm(max.tvtr)
   rm(min.tvtr)
   
   cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
   cat("\nTotal number of samples to exclude because of TR-TV ratio is",length(cols.to.exclude),"\n.....")
   #log file
   cat(paste("Total number of samples to exclude because of TR-TV ratio is: ",length(cols.to.exclude)),
       file=logFile, append=TRUE, sep = "\n") 
   
   if (length(cols.to.exclude) > 0){
      mut.all <- mut.all[,-cols.to.exclude]   
      melt.mut.all <- melt.mut.all[!melt.mut.all$variable %in% samples.to.exclude,]
      melt.mut.all <- droplevels(melt.mut.all)
      samples.info2 <- samples.info2[!samples.info2[,answ.name] %in% samples.to.exclude,]
      samples.info2 <- droplevels(samples.info2)
      #all.samples <- all.samples[!all.samples %in% samples.to.exclude]
      #recaluclate annotation and sample columns
      annot.cols <- which(grepl('#',colnames(mut.all)) )
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      
      remove.snps <- which(rowSums(mut.all[,samples.cols], na.rm=T) == 0)
      if(length(remove.snps)>0){
         mut.all <- mut.all[-remove.snps,]  
         cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
         cat(paste("Mutations with 0 (or only NAs) samples mutated are excluded. Amount of such mutations is:",length(remove.snps), sep=' '), 
             file=logFile, append=TRUE, sep = "\n") }
      rm(remove.snps) 
   }
   rm(cols.to.exclude)
 
} else {
   rm(samples.to.exclude.temp)
   rm(max.tvtr)
   rm(min.tvtr)
   cat("You choose not to remove any sample. \n.....")
}

melt.mut.all$tvtr <- NULL

#log file
cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of samples left after filtering:",nrow(samples.info2)),  file=logFile, append=TRUE, sep = "\n")
#save.image('Temp.data4.RData')

#################################################################################
# to plot in the end samples number of mutations distribution (colored by type)
#################################################################################

if (exists('col.snps.nar')){
   cat("\n ---------------------------------------------------------------------------------------------\n")
   cat("\n --------- Chooseing grouping and naming for samples number of mutations in barplot ----------\n")
   cat("\n ---------------------------------------------------------------------------------------------\n")
   
   # lof file
   cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
   cat('---------- Chooseing grouping and naming for samples number of mutations in barplot ----------', file=logFile, append=TRUE, sep = "\n")
   cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")
   
   #log file
   cat(paste("Barplot for number of mutations per samples will use", col.snps.nar,"to color group of mutation types"),  file=logFile, append=TRUE, sep = "\n")
   cat("------------------------------------------",  file=logFile, append=TRUE, sep = "\n")
   
   options <- as.character(unique(melt.mut.all[,col.snps.nar]))
   options[is.na(options)] <- 'NA_Missing_info'
   cat("\n\nIn column",col.snps.nar,"there are folowing fields:")
   counts.tb <- table(melt.mut.all[,col.snps.nar], useNA='always')
   names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
   counts.tb <- counts.tb[options]
   cat("\n",paste(1:length(options),".", options ," (", counts.tb ," mutations)","\n", sep=''))
   cat("Do you want to merge and  give new name to some of the columns, for coloring purposes? Answer y or n : ")
   col.ind <- scan(read.from,what=character(),n=1,quiet=TRUE, nlines=1)
   col.ind <- tolower(substr(col.ind,1,1))
   while(col.ind=='y'){
      cat("Which one you want to combine: ")
      combine.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
      combine.options <- unique(combine.options)
      combine.options <- combine.options[combine.options <= length(options)]
      cat("How you would like to name them: ")
      rename.options <- scan(read.from,what=character(),nmax=1,quiet=TRUE, nlines=1)
      if( 'NA_Missing_info' %in% options[combine.options]){
         melt.mut.all[is.na(melt.mut.all[,col.snps.nar]),col.snps.nar] <- rename.options   
         if (length(combine.options) > 1){
            melt.mut.all[melt.mut.all[,col.snps.nar] %in% options[combine.options], col.snps.nar]<- rename.options
         }
      } else {
         melt.mut.all[melt.mut.all[,col.snps.nar] %in% options[combine.options] &
                       !is.na(melt.mut.all[,col.snps.nar]), col.snps.nar]<- rename.options      
      }
      cat(paste("Mutations which have", col.snps.nar,"matching:",paste0(options[combine.options],collapse=', '),
                " - will be renamed to",rename.options), 
          file=logFile, append=TRUE, sep = "\n")
      
      ##
      options <- as.character(unique(melt.mut.all[,col.snps.nar]))
      options[is.na(options)] <- 'NA_Missing_info'
      if(length(options)>1){
         cat("\n\nIn column",col.snps.nar,"there are folowing fields now:")
         counts.tb <- table(melt.mut.all[,col.snps.nar], useNA='always')
         names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
         counts.tb <- counts.tb[options]
         cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' mutations)','\n', sep=''))
         cat("Do you want to merge more? Answer y or n : ")
         col.ind <- scan(read.from,what=character(),n=1,quiet=TRUE, nlines=1)
         col.ind <- tolower(substr(col.ind,1,1))
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
#rm(char.cols)

ordered.samples <- as.character(names(sort(table(melt.mut.all$variable), decreasing=T)))
melt.mut.all <- melt.mut.all[,c('variable', col.snps.nar)]
colnames(melt.mut.all) <- c('Sample', 'Type')
temp.df <- table(melt.mut.all[,c('Sample','Type')], useNA="always")
melt.fin.df <- melt(temp.df)

melt.fin.df$Sample <- factor(melt.fin.df$Sample, levels=ordered.samples)

rm(temp.df)
rm(melt.mut.all)

samples.info2[,answ.name] <- factor(samples.info2[,answ.name], levels=ordered.samples)
#save.image('Temp.data5.RData')

######################################################################################
# PCA QC and filtering
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n ----------------- Filtering outlier samples based on PCA 1 and 2 components -----------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('------------------ Filtering outlier samples based on PCA 1 and 2 components -----------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))
snpid.col <- match(c('#snpID','#Chr','#Position','#Reference','#Alteration'), colnames(mut.all))

cat("\n\nYou  want to do PCA on filtered mutations? \nChoose y if yes, or n to do it on all mutations(y or n): ")
answ.pca <- scan(read.from,what=character(),n=1,quiet=TRUE)
answ.pca <- tolower(substr(answ.pca,1,1))
if(answ.pca=='y'){    
   cols.types <- sapply(mut.all[,annot.cols], class)
   char.cols <- cols.types[cols.types =='character']
   rm(cols.types)
   cat("\n\nWhich column you want to use to filter:", 
       paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
   cat("\nEnter colum number: ")
   x.col <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   x.col <- unique(x.col)
   x.col <- x.col[x.col<= length(char.cols)] 
   if (length(x.col) > 0){
      filt.pca.snps <- names(char.cols[x.col])
      cat("\nColumn you chose is:",filt.pca.snps, "\n.....")
      cat(paste("Mutations used for PCA are filterd by:", filt.pca.snps),  file=logFile, append=TRUE, sep = "\n")
      options <- as.character(unique(mut.all[,filt.pca.snps]))
      options[is.na(options)] <- 'NA_Missing_info'
      cat("\n\nIn column",filt.pca.snps,"there are folowing fields:")
      counts.tb <- table(mut.all[,filt.pca.snps], useNA='always')
      names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
      counts.tb <- counts.tb[options]
      cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' mutations)','\n', sep=''))
      cat("Whish one you wnat to use for PCA plot: ")
      pca.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
      pca.options <- unique(pca.options)
      pca.options <- pca.options[pca.options <= length(options)]
      if( 'NA_Missing_info' %in% options[pca.options]){
         filter.mut <- mut.all[ mut.all[,filt.pca.snps ] %in%  options[pca.options] |  is.na(mut.all[,filt.pca.snps]), c(snpid.col, samples.cols)]          
      } else {
         filter.mut <- mut.all[ mut.all[,filt.pca.snps] %in% options[pca.options] &  !is.na(mut.all[,filt.pca.snps]), c(snpid.col,samples.cols)]                
      }
      cat(paste("Mutations used for PCA are one that have", filt.pca.snps, "matching:",paste0(options[pca.options], collapse=', ' )),
          file=logFile, append=TRUE, sep = "\n")
      
      rm(pca.options)
      rm(counts.tb)
      rm(options)
      rm(filt.pca.snps)
   } else {
      cat("\nYou choose not to filter in the end.\n.....")
      cat(paste("No filtering is performed prior to PCA calculation."),  file=logFile, append=TRUE, sep = "\n")
      filter.mut <- mut.all[, c(snpid.col,samples.cols)]
   }
   rm(char.cols)
   rm(x.col)    
} else {
   cat("\nPCA will be done on all (not linked) mutations.\n.....")
   cat(paste("No filtering is performed prior to PCA calculation."),  file=logFile, append=TRUE, sep = "\n")
   #not realy silent, just keep name for siplicity
   filter.mut <- mut.all[,c(snpid.col, samples.cols)] 
} 
snpid.col <-  which(grepl('#',colnames(filter.mut)))

#only rare snps option
cat("\nYou have",ncol(filter.mut)-length(snpid.col),"patients, and ",nrow(filter.mut),"snps for PCA. You want to remove rare snps (<0.005) ? (y or n): ")
answ.pca.rare <- scan(read.from,what=character(),n=1,quiet=TRUE)
answ.pca.rare <- tolower(substr(answ.pca.rare,1,1))
if (answ.pca.rare == 'y'){
   # remove singltons
   # mat.pat <- mat.pat[,colSums(mat.pat, na.rm=T)>1]  
   a <- as.numeric(rep(NA, nrow(filter.mut)))
   a <- apply(filter.mut[,-snpid.col],1,function(x) sum(x, na.rm=T)/(sum(!is.na(x))*2))
   filter.mut <- filter.mut[a>0.005,]
   rm(a)
   cat(paste("Only mutations with AF above 0.005 are kept for PCA."),  file=logFile, append=TRUE, sep = "\n")
   cat("\nOnly mutations with AF above 0.005 are kept for PCA.\n.....")  
}
rm(answ.pca.rare)
#save.image('Temp.data6.RData')

#get only sample cols
#pca.samples.cols <- which(!grepl('#',colnames(filter.mut)))
# exclude samples with no patient mutated
#filter.mut <- filter.mut[ ,colSums(filter.mut[,pca.samples.cols], na.rm=T) != 0]
#rm(pca.samples.cols)

#snps info
snp.info <-filter.mut[,snpid.col]
# matrix with genotypes only (0,1,2), no annotation columns
filter.mut <- (as.matrix(filter.mut[,-(snpid.col)]))

rm(snpid.col)

###################################################
## Pruning using SNPRelate
###################################################
sample.id <- colnames(filter.mut)
snp.id <- snp.info$'#snpID'
snp.chr <- snp.info$'#Chr'
snp.pos <- snp.info$'#Position'
snp.allele <- paste(snp.info$'#Reference',snp.info$'#Alteration',sep='/')

snpgdsCreateGeno("casecon.gds", genmat=filter.mut, sample.id=sample.id, snp.id=snp.id,
                 snp.chromosome=snp.chr, snp.position=snp.pos, snp.allele=snp.allele, snpfirstdim=TRUE)

gdsobj <- snpgdsOpen("casecon.gds")

snpset <- snpgdsLDpruning(gdsobj,ld.threshold = 0.2)
snp.sel <- unlist(snpset)
snp.sel.num <- which(snp.id %in% snp.sel)
mat.pat <- t(filter.mut[ snp.sel.num , ])

cat(paste("Number of mutations before LD filter, for PCA calculation:", length(snp.id)),  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of mutations with no LD used for PCA calculation:", ncol(mat.pat)),  file=logFile, append=TRUE, sep = "\n")

snpgdsClose(gdsobj)
rm(snp.id)
rm(snp.chr)
rm(snp.pos)
rm(snp.allele)
rm(snp.info)
rm(snp.sel)
rm(snp.sel.num)
rm(sample.id)
rm(filter.mut)
rm(snpset)
rm(gdsobj)

################################

# reduce if too much
if(ncol(mat.pat)>100000){
   mat.pat <- mat.pat[,sample(1:ncol(mat.pat), size=100000, replace=F)]
   cat("\nPCA will be done on random 100000 mutaations (after pruning).\n.....")
   cat(paste("PCA will be done on random 100000 mutaations (after pruning)"),  file=logFile, append=TRUE, sep = "\n")
} else {
   cat("\nPCA will be done on",ncol(mat.pat),"mutations (after pruning). And",nrow(mat.pat),"samples.\n.....")
   cat(paste("PCA will be done on",ncol(mat.pat)," mutations (after pruning). And",nrow(mat.pat),"samples."),  file=logFile, append=TRUE, sep = "\n")  
}

# choose color and shape columns
cat("\n\nYou can color/shape in PCA plot by:", 
    paste('\n',1:(ncol(samples.info2)),'.',colnames(samples.info2), sep=''))
cat("\nEnter  numbers space separated for color and shape (in this order, or one/none if you don't want color/shape): ")
x <- scan(read.from,what=integer(),nmax=ncol(2),quiet=TRUE, nlines=1)
if (length(x) > 0){
   color.col <- colnames(samples.info2)[x[1]]
   cat(paste("Samples color at PCA plot is by:",color.col),  file=logFile, append=TRUE, sep = "\n")  
} 

if (length(x) > 1){
   shape.col <- colnames(samples.info2)[x[2]]
   cat(paste("Samples shape at PCA plot is by:",shape.col),  file=logFile, append=TRUE, sep = "\n")  
} 

if (exists('shape.col')){
   num.shapes <- length(unique(samples.info2[,x[2]]))
   if (num.shapes > 25 ){
      vv <- c(sample(1:25),sample(1:25,size=(num.shapes-25), replace=F))   
   } else {
      vv <- c(sample(1:num.shapes))
   }
   rm(num.shapes)
}
#rm(x)

pc.max <- min(20,min(dim(mat.pat)))

#replace NAs for variance estimatece
if(sum(is.na(mat.pat)) >0 ){
   mat.pat2 <- mat.pat
   mat.pat2[is.na(mat.pat2)] <- 0
   temp.pca <- pca(X=mat.pat2, ncomp=pc.max)
   rm(mat.pat2)
} else {
   temp.pca <- pca(X=mat.pat, ncomp=pc.max)
}

#var.nip <- tune.pca(X=mat.pat2, center=F, max.iter=300 , ncomp=min(20,min(dim(mat.pat2))), tol=1e-07)
#dev.off()


dd <- data.frame(Var=temp.pca$explained_variance, Cum.Var=temp.pca$cum.var, comp=names(temp.pca$explained_variance))
dd$comp <- factor(dd$comp, levels=names(temp.pca$explained_variance))
pos.samp <- which(answ.name==colnames(samples.info2))
ddr <- merge(temp.pca$x, samples.info2[,c(pos.samp, x)], by.y=answ.name, by.x="row.names")

if (length(x) > 0){
   ddr[,color.col] <- factor(ddr[,color.col])
} else {
   color.col <- 'color_dummy'
   ddr$color_dummy <- factor('None')
}

p.1 <- ggplot(dd, aes(y=Var, x=comp)) + geom_bar(stat = "identity") + theme_bw() +
   ylab(" % Exp. Var. (aprox)") + xlab("Components")
p.2 <- ggplot(dd, aes(y=Cum.Var, x=comp)) + geom_point() + geom_line(group=1) + theme_bw() + 
   ylab("Cumulative % EV (aprox)") + xlab("Components")
p.3 <- ggpairs(ddr[,-1], 1:min(10,pc.max), mapping=ggplot2::aes_string(colour=color.col, alpha=0.5), legend = c(2,1),
               lower = list(continuous = wrap("points", alpha = 0.75, size=1)), 
               upper = list(continuous = 'blank')) +
               theme_bw() + theme(legend.position = 'right') 
g.3 <- grid.grabExpr(print(p.3))
X11(width=14, height=8)
#par(mfrow=c(2,1))
grid.arrange(p.1, p.2, g.3, layout_matrix = rbind(c(1,1,2,2),c(3,3,3,3),c(3,3,3,3),c(3,3,3,3)))



#plot(var.nip$cum.var, type='l', ylim=c(0,1), col=2, xlab="PC components", ylab="Proportion of Explaine Variance")
#barplot(var.nip$explained_variance, ylab="% of Explained Variance (aproximation)")
#plot(var.nip$cum.var, type='o', pch=19, xlim=c(0,pc.max), 
#     ylim=c(min(var.nip$cum.var[1:min(20,pc.max)])*0.98, max(var.nip$cum.var[1:min(20,pc.max)])*1.02),
#     xlab="PC components (first 20)", ylab="Cumulative % of Exp. Var. (aproximation)")

# ask user how many components tey want
cat("\nEnter  numbers of principal components you want to be calculated (n <= 20) \n(Warrning: only first two PC comp. will be plotted for filtering): ")
n.cp <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
n.cp <- min(c(n.cp, 20))
dev.off()

pdf(file='PCA_variance_explained_old.pdf', width=14,height=14,useDingbats=F)
grid.arrange(p.1, p.2, g.3, layout_matrix = rbind(c(1,1,2,2),c(3,3,3,3),c(3,3,3,3),c(3,3,3,3)))
dev.off()
#log file
cat(paste("PCA_variance_explained_old.pdf done! - PCA first two componets - not final"),
    file=logFile, append=TRUE, sep = "\n")


rm(temp.pca)
rm(p.1)
rm(p.2)
rm(p.3)
rm(g.3)
#rm(var.nip)
cat(paste("Number of PCA components:", n.cp),  file=logFile, append=TRUE, sep = "\n")

# first 2 copmponents of PCA
mat.pat.nip <- pca(X=mat.pat, ncomp=2)
df <- data.frame( sample_ID = rownames(mat.pat.nip$x),
                  PC1 = mat.pat.nip$x[,1],
                  PC2 = mat.pat.nip$x[,2])

rm(mat.pat.nip)

pca.snps.all <- merge(x=df, y=samples.info2, by.x='sample_ID', by.y=answ.name,  all.x=T)
rm(df)
#rm(mat.pat)

######

cat("\n\nIn the plot choose cordinates in which you want to keep the data. \n
    (is there data points which look like outliers? Choose thresholds for x-min, x-max, y-min and y-max) \n .....")


X11(width=14, height=8)
if (exists('color.col') & exists('shape.col')){
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', colour=color.col, shape=shape.col )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
   p2
} else if (exists('color.col')) {
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', colour=color.col )) +
      geom_point(size=4, alpha=0.8) + theme_bw()
   p2
} else if (exists('shape.col')) {
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', shape=shape.col )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
   p2
} else {
   p2 <-ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2' )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
   p2
}

cat("\nEnter x-min  value (just press enter if you don't want any): ")
xmin <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter x-max value(just press enter if you don't want any): ")
xmax <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter y-min value(just press enter if you don't want any): ")
ymin <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter y-max value(just press enter if you don't want any): ")
ymax <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
#log file
cat(paste("Minimum PC1 component:",xmin),  file=logFile, append=TRUE, sep = "\n")
cat(paste("Maximum PC1 component:",xmax),  file=logFile, append=TRUE, sep = "\n")
cat(paste("Minimum PC2 component:",ymin),  file=logFile, append=TRUE, sep = "\n")
cat(paste("Maximum PC2 component:",ymax),  file=logFile, append=TRUE, sep = "\n")

dev.off()

pdf(file='projects.PC1.PC2.pca_old.pdf', width=20,height=20,useDingbats=F)
p2
dev.off()
rm(p2)
#log file
cat(paste("projects.PC1.PC2.pca_old.pdf done! - PCA first two componets - not final"),
    file=logFile, append=TRUE, sep = "\n")

samples.to.exclude.temp <- c()
if (length(xmin) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC1 < xmin,'sample_ID']))
}
if (length(xmax) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC1 > xmax,'sample_ID']))
}
if (length(ymin) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC2 < ymin,'sample_ID']))
}
if (length(ymax) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC2 > ymax,'sample_ID']))
}
samples.to.exclude.temp <- unique(samples.to.exclude.temp)

if(length(samples.to.exclude.temp) >0){
   cat("Samples (",length(samples.to.exclude.temp)," of them ) you decided to exclude because of PCA plot are:",paste(samples.to.exclude.temp, collapse=", "), "\n .....")
   #log file
   cat(paste("Samples exclude because of PCA are:", paste(samples.to.exclude.temp, collapse=", ")),
       file=logFile, append=TRUE, sep = "\n") 
   
   samples.to.exclude <- c(samples.to.exclude, samples.to.exclude.temp)
   rm(samples.to.exclude.temp)
   
   cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
   cat("\nTotal number of samples to exclude because of PCA is",length(cols.to.exclude), "\n .....")
   #log file
   cat(paste("Total number of samples to exclude because of PCA is: ",length(cols.to.exclude)),
       file=logFile, append=TRUE, sep = "\n") 
   
   if (length(cols.to.exclude) > 0){
      mut.all <- mut.all[,-cols.to.exclude]
      mat.pat <- mat.pat[!rownames(mat.pat) %in% samples.to.exclude,]
      mat.pat <- mat.pat[,colSums(mat.pat, na.rm=T)!=0]
       
      melt.fin.df <- melt.fin.df[!melt.fin.df$Sample %in% samples.to.exclude,]
      melt.fin.df$Sample <- droplevels(melt.fin.df$Sample)
      
      samples.info2 <- samples.info2[!samples.info2[,answ.name] %in% samples.to.exclude,]
      samples.info2 <- droplevels(samples.info2)
      
      #ordered.samples <- ordered.samples[!ordered.samples %in% samples.to.exclude ]
      #samples.info2[,answ.name] <- factor(samples.info2[,answ.name], levels=ordered.samples)
      
      annot.cols <- which(grepl('#',colnames(mut.all)) )
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      
      remove.snps <- which(rowSums(mut.all[,samples.cols], na.rm=T) == 0)
      if(length(remove.snps)>0){
         mut.all <- mut.all[-remove.snps,]  
         cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
         cat(paste("Mutations with 0 samples mutated are excluded. Amount of such mutations is:",length(remove.snps), sep=' '), 
             file=logFile, append=TRUE, sep = "\n") 
      }
      rm(remove.snps) 
      
   }
   rm(cols.to.exclude) 
   
} else {
   rm(samples.to.exclude.temp)  
   cat("\nNo samples are chosen to be remove because of PCA outliers\n .....")
}
rm(pca.snps.all)
rm(xmin)
rm(xmax)
rm(ymax)
rm(ymin)
rm(pc.max)
rm(dd)
rm(ddr)
#log file
cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of samples left after filtering:",nrow(samples.info2)),  file=logFile, append=TRUE, sep = "\n")
#save.image('Temp.data7.RData')


######################################################################################
# Repolot everything after final filtration
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n --------------- Recalulating, Reploting, and saving output for further steps ----------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('---------------- Recalulating, Reploting, and saving output for further steps ----------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

if(exists('col.snps.nar')){
   p1 <- ggplot(data=melt.fin.df, aes(x=Sample, y=value, fill=Type)) + geom_bar(stat="identity") + 
      theme_bw() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
      scale_y_continuous(expand=c(0,0)) + xlab('Samples') + ylab('Number of mutations')+ 
      ggtitle('Barplot of number of mutations per sample, colored by type ')
   pdf(file='sample.mutation.distribution.pdf', width=25,height=15,useDingbats=F)
   print(p1)
   dev.off()  
   #log file
   cat(paste("sample.mutation.distribution.pdf done! - Barplot by type of mutations - final"),
       file=logFile, append=TRUE, sep = "\n")
   rm(p1)
   rm(melt.fin.df)
   rm(col.snps.nar) 
}

#plot hist of per saple number of mutations
pdf(file='hist.number.mutations.pdf', width=12,height=12,useDingbats=F)
ggplot(samples.info2, aes(x=Num.Mutations)) + geom_histogram() + theme_bw() + ggtitle('Number of mutations per sample histogram')
dev.off()
#log file
cat(paste("hist.number.mutations.pdf done! - histogram for TR-TV ratio per sample - final"),
    file=logFile, append=TRUE, sep = "\n")

pdf(file='hist.trtv.ratio.pdf', width=12,height=12,useDingbats=F)
ggplot(samples.info2, aes(x=tvtr_ratio)) + geom_histogram() + theme_bw()  + ggtitle('Tr/Tv ratio per sample histogram')
dev.off()
#log file
cat(paste("hist.trtv.ratio.pdf done! - Barplot by type of mutations - final"),
    file=logFile, append=TRUE, sep = "\n")
##############
# repeat PCA
##############
# first n copmponents of PCA
n.cp2 <- max(2,n.cp)
pca.fin <- pca(X=mat.pat, ncomp=n.cp2)

dd <- data.frame(Var=pca.fin$explained_variance, Cum.Var=pca.fin$cum.var, comp=names(pca.fin$explained_variance))
dd$comp <- factor(dd$comp, levels=names(pca.fin$explained_variance))
pos.samp <- which(answ.name==colnames(samples.info2))
pca.snps.all <- merge(pca.fin$x,  samples.info2, by.y=answ.name, by.x="row.names")

p.1 <- ggplot(dd, aes(y=Var, x=comp)) + geom_bar(stat = "identity") + theme_bw() +
   ylab(" % Exp. Var. ") + xlab("Components")
p.2 <- ggplot(dd, aes(y=Cum.Var, x=comp)) + geom_point() + geom_line(group=1) + theme_bw() + 
   ylab("Cumulative % Exp. Var.") + xlab("Components")
if (length(x) > 0){
   p.3 <- ggpairs(pca.snps.all[,-1], 1:n.cp2, mapping=ggplot2::aes_string(colour=color.col, alpha=0.5), legend = c(2,1),
                  lower = list(continuous = wrap("points", alpha = 0.75, size=1)), 
                  upper = list(continuous = 'blank')) +
      theme_bw() + theme(legend.position = 'right') 
} else {
   p.3 <- ggpairs(pca.snps.all[,-1], 1:n.cp2, mapping=ggplot2::aes(alpha=0.5), 
                  lower = list(continuous = wrap("points", alpha = 0.75, size=1)), 
                  upper = list(continuous = 'blank')) +  theme_bw()
}
g.3 <- grid.grabExpr(print(p.3))

pdf(file='PCA_variance_explained.pdf', width=14,height=14,useDingbats=F)
grid.arrange(p.1, p.2, g.3, layout_matrix = rbind(c(1,1,2,2),c(3,3,3,3),c(3,3,3,3),c(3,3,3,3)))
dev.off()

rm(pca.fin)

#log file
cat(paste("PCA_variance_explained.pdf done! - Prcenteges of variance explained by PCA components - final"),
    file=logFile, append=TRUE, sep = "\n")
rm(p.1)
rm(p.2)
rm(p.3)
rm(g.3)

#plot pca
if (exists('color.col') & exists('shape.col')){
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', colour=color.col, shape=shape.col )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
} else if (exists('color.col')) {
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', colour=color.col )) +
      geom_point(size=4, alpha=0.8) + theme_bw()
} else if (exists('shape.col')) {
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', shape=shape.col )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
} else {
   p2 <-ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2' )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
}
pdf(file='projects.PC1.PC2.pca.pdf', width=18,height=15,useDingbats=F)
p2
dev.off()
rm(p2)
#log file
cat(paste("projects.pca.pdf done! - PCA first two componets - final"),
    file=logFile, append=TRUE, sep = "\n")

rm(mat.pat)
#rm(mat.pat.nip)
#rm(pc.max)

cat("\nThere is left ", nrow(mut.all)," mutations and indels together ")
#log file
cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of mutations left:", nrow(mut.all)),file=logFile, append=TRUE, sep = "\n")

cat("\nTable dimensions are: ", dim(mut.all), "(",ncol(mut.all)-length(annot.cols)," samples ) \n.....")
#log file
cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
cat(paste("Number of samples left:",nrow(pca.snps.all)),file=logFile, append=TRUE, sep = "\n")

write.table(mut.all, file='mutations_filtered_samples.txt', quote=F, sep='\t', row.names=F)
cat("\n *****************************************************************************************************************\n")
cat("\n * File you should use in next steps of analysis as mulisample mutations file is: mutations_filtered_samples.txt *")
cat("\n *****************************************************************************************************************\n")
#log file
cat("\n******************************************************************************************",file=logFile, append=TRUE, sep = "\n")
cat(paste("*        mutations_filtered_samples.txt file done! Use this file in further steps        *"), file=logFile, append=TRUE, sep = "\n")
cat("******************************************************************************************",file=logFile, append=TRUE, sep = "\n")
#save(mut.all, file='mut.all.RData')
# correct it only one PCA is chosen
if(n.cp2 != n.cp) {
   df <- df[,1:(n.cp+1)] 
   pca.snps.all <- merge(x=df, y=samples.info2,  by.x='sample_ID', by.y=answ.name, all.x=T)
}

cols.types <- sapply(pca.snps.all[,-1], class)
num.cols <- cols.types[cols.types =='numeric' |cols.types =='integer' ]
num.cols <- names(num.cols)
not.num <- cols.types[cols.types !='numeric' & cols.types !='integer' ]
not.num <- names(not.num)
sample.id <- colnames(pca.snps.all)[1]
pca.snps.all <- pca.snps.all[,c(sample.id, num.cols, not.num )]
rm(cols.types)
rm(num.cols)
rm(not.num)
rm(sample.id)

write.table(pca.snps.all, file='samples_info_pca.txt', quote=F, sep='\t', row.names=F)
cat("\n *****************************************************************************************************************\n")
cat("\n *   File you should use in next steps of analysis as covariate file (samples info file): samples_info_pca.txt   *")
cat("\n *****************************************************************************************************************\n")
#log file
cat("\n******************************************************************************************",file=logFile, append=TRUE, sep = "\n")
cat(paste("*    samples_info_pca.txt file done! Use this file for covariates in further steps       *"), file=logFile, append=TRUE, sep = "\n")
cat("******************************************************************************************",file=logFile, append=TRUE, sep = "\n")

rm(n.cp)
rm(n.cp2)
rm(vv)
rm(shape.col)
rm(color.col)
rm(answ.pca)
save.image('Temp.data7.RData')
######################################################################################
# distribution cases vs controls
######################################################################################
cat("\n ---------------------------------------------------------------------------------------------\n")
cat("\n --------------------- (Optional) Case vs. Contorls explorative analysis ---------------------\n")
cat("\n ---------------------------------------------------------------------------------------------\n")

# lof file
cat('\n----------------------------------------------------------------------------------------------', file=logFile, append=TRUE, sep = "\n")
cat('---------------------- (Optional) Case VS Contorls explorative analysis ----------------------', file=logFile, append=TRUE, sep = "\n")
cat('----------------------------------------------------------------------------------------------\n', file=logFile, append=TRUE, sep = "\n")

## compare mutation distribution of two different cohorts
cat("\n\nYou  want to compare two cohorts ? you can split them by one of the following columns", 
    paste('\n',1:ncol(samples.info2),'.',colnames(samples.info2), sep=''),
    "\nYou want to comapare (y or n): ")
answ <- scan(read.from,what=character(),n=1,quiet=TRUE)
if(answ=='y'){  
   cat("\nEnter column number of column you want ot use to split to controls and cases: ")
   x <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   x <- x[x<= ncol(samples.info2)]       
   
   cat("\nColumn(s) you chose is:",paste(colnames(samples.info2)[x], collapse=', ' ), "\n.....")   
   options <- as.character(unique(samples.info2[,x]))
   options[is.na(options)] <- 'NA_Missing_info'
   cat("\nIn column",colnames(samples.info2)[x],"there are folowing fields:")
   counts.tb <- table(samples.info2[,x], useNA='always')
   names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
   counts.tb <- counts.tb[options]
   cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' samples)','\n', sep=''))
   cat("\nEnter  numbers of options you consider as cases (space separated): ")
   cases.name <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
   cases.name <- unique(cases.name)
   cases.name <- cases.name[cases.name <= length(options)]
   
   if (length(cases.name) > 0) {
      cat("\nSamples whos", colnames(samples.info2)[x],
          "is in:", paste(options[cases.name], collapse=', '),
          "- will be considered as cases")      
      # log file
      cat(paste("Samples that have :",colnames(samples.info2)[x], "matching:",paste(options[cases.name], collapse=', '),
                "- will be considered as cases"),file=logFile, append=TRUE, sep = "\n")
      cases <-  as.character(samples.info2[as.character(samples.info2[,x]) %in% options[cases.name],answ.name])
      cases <- unique(cases) 
   } else {
      stop(" You didn't choose any option, you should have said no (n) before...\n ")
   }     
   cat("\n ---------------------------------------------------------------------------------------------\n")
   rm(cases.name)         
   rm(options)
   rm(counts.tb)
   rm(x)
   
   samples.info2$status <- 'controls' 
   samples.info2[ samples.info2[,answ.name] %in% cases ,'status'] <- 'cases'
   
   samples.info2[,answ.name] <- factor(samples.info2[,answ.name], 
                                       levels=unique(as.character(samples.info2[order(-samples.info2$Num.Mutations) ,answ.name])))
   
   p1 <- ggplot(data=samples.info2, aes_string(x=answ.name, y='Num.Mutations', fill='status')) + geom_bar(stat="identity")   + 
      theme_bw() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
      scale_y_continuous(expand=c(0,0)) + xlab('Samples') + ylab('Number of mutations')
   #save(p1, file='plot_sample.mutation.hist_caseVScontrol.RData')
   
   pdf(file='samples.hist_caseVScontrol.pdf', width=25,height=15,useDingbats=F)
   print(p1)
   dev.off() 
   #log file
   cat(paste("samples.hist_caseVScontrol.pdfe! - Barplot for number of mutations per sample, colored by cases and contorls - final"),
      file=logFile, append=TRUE, sep = "\n")

 
   rm(p1)
}



if(answ=='y' & TRUE ){
   #density of mutations per gene
   cases.num <- sum(samples.info2$status == 'cases')
   #cases.samp <- as.character(samples.info2[samples.info2$status == 'cases','samples'])
   cases.samp <- cases
   control.num <-  sum(samples.info2$status == 'controls')
   controls.samp <- as.character(samples.info2[samples.info2$status == 'controls',answ.name])
   
   # log file
   cat(paste("There are",cases.num,"cases"),file=logFile, append=TRUE, sep = "\n")
   cat(paste("There are",control.num,"controls"),file=logFile, append=TRUE, sep = "\n")

   # choose column to aggregate on , e.g. Genes
   annot.cols <- which(grepl('#',colnames(mut.all)))
   samples.cols <- which(!grepl('#',colnames(mut.all)))

   # Aggregate column
   cols.types <- sapply(mut.all[,annot.cols], class)
   num.cols <- cols.types[cols.types =='integer' | cols.types =='numeric']
   char.cols <- cols.types[cols.types =='character']
   rm(cols.types)
   cat("\n\nYou can choose chaarachter  colums to aggregate on (e.g. Genes):", 
       paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
   cat("\nEnter column number: ")
   x <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
   x <- unique(x)
   x <- x[x<= length(char.cols)]    
   if (length(x) <= 0){
      stop("\nYou  choose wrong columns .....")  
   } else {
      agg.name <- names(char.cols)[x]
      cat("\nColumn(s) you chose is:",agg.name,  "\n.....")
      # log file
      cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
      cat(paste("Aggregation of mutations with be done on",agg.name,"column"),file=logFile, append=TRUE, sep = "\n")
   }  
   rm(x)
  
   cat("\n ---------------------------------------------------------------------------------------------\n")
   #log file
   cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")

   # filtering chategorical column
   cat("\n\nYou can choose chaarachter  colums to filter by (e.g Exonic function, only silent snps):", 
       paste('\n',1:length(char.cols),'.',names(char.cols), sep='' ))
   cat("\nEnter the column numbers, space separated: ")
   x <- scan(read.from,what=integer(),nmax=length(char.cols),quiet=TRUE, nlines=1)
   x <- unique(x)
   x <- x[x<= length(char.cols)]       
   if (length(x) <= 0){
      cat("\nYou or choose not to filter or wrong columns \n.....")
      cat(paste("For case-control compariosn non of mutations will be filtered besed on character annotation columns"),
          file=logFile, append=TRUE, sep = "\n")
   } else {
      cat("\nColumn(s) you chose is:",paste(names(char.cols)[x], collapse=', ' ), "\n.....")
      for (i in x){
         temp.col <- names(char.cols[i])
         options <- as.character(unique(mut.all[,temp.col]))
         options[is.na(options)] <- 'NA_Missing_info'
         cat("\n\nIn column",names(char.cols)[i],"there are folowing fields:")
         counts.tb <- table(mut.all[,temp.col], useNA='always')
         names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
         counts.tb <- counts.tb[options]
         cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' mutations)','\n', sep=''))
         cat("Enter the numbers which are in front of the options you DON'T want to keep (space separated): ")
         exclude.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
         exclude.options <- unique(exclude.options)
         exclude.options <- exclude.options[exclude.options <= length(options)]
         if (length(exclude.options) > 0) {
            cat("\nMutations whos", temp.col,
                "is in:", paste(options[exclude.options], collapse=', '),
                "- will be excluded")  
            # log file
            cat(paste("Mutations that have",temp.col,"matching",paste(options[-exclude.options], collapse=', '),
                      "- will be only keept"),file=logFile, append=TRUE, sep = "\n")
            if( 'NA_Missing_info' %in% options[exclude.options]){
               mut.all <-  mut.all[!is.na(as.character(mut.all[,temp.col])) , ]
            }
            mut.all <- mut.all[! mut.all[,temp.col] %in% options[exclude.options] | is.na(mut.all[,temp.col]),  ]
         }  else {
            cat("\n Nothing from column",temp.col,"will be excluded ")
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
   # filtering numeric column
   cat("\n\nYou can choose numeric colums to filter by:", 
       paste('\n',1:length(num.cols),'.',names(num.cols), sep='' ))
   cat("\nEnter the column numbers, space separated: ")
   x <- scan(read.from,what=integer(),nmax=length(num.cols),quiet=TRUE, nlines=1)
   x <- unique(x)
   x <- x[x<= length(num.cols)]       
   if (length(x) <= 0){
      cat("\nYou or chose not to filter or wrong columns \n.....")
      cat(paste("For case-control compariosn non of mutations will be filtered besed on numeric annotation columns"),
          file=logFile, append=TRUE, sep = "\n")      
   } else {
      cat("\nColumn(s) you chose is:",paste(names(num.cols)[x], collapse=', ' ), "\n.....")
      for (i in x){
         temp.col <- names(num.cols[i])
         quant.temp.col <- quantile(mut.all[,temp.col], na.rm=T)
         NA.temp.col <- sum(is.na(mut.all[,temp.col]))
         if(NA.temp.col >0 ){
            cat("\nIn column",temp.col,"therea are",NA.temp.col,'Na values')
            cat("\n\nDo you want filter out these mutations (answer y or n)? \nAnswer:")
            NA.ind <- scan(read.from,what=character(),n=1,quiet=TRUE)
            if (NA.ind=='y'){
               mut.all <- mut.all[!is.na(mut.all[,temp.col]),  ]
            }
            rm(NA.ind)                   
         }
         cat("\nIn column",temp.col,"quantiles values are:", 
             paste(paste0('\n',quant.temp.col), names(quant.temp.col), sep=' at qunatile ' ) )
         cat("\n\nDo you want filter higher or lower then treshold (h for higer, l for lower)? \nAnswer: ")
         num.ind <- scan(read.from,what=character(),n=1,quiet=TRUE)
         num.ind <- tolower(num.ind)
         num.ind <- substr(num.ind,1,1)
         cat("\nWhat is treshold you want to use: ")
         trehsold <- scan(read.from,what=numeric(),n=1,quiet=TRUE)       
         if(num.ind == 'l'){
            message('Everything lower then ',trehsold,' in column ', temp.col, ' will be removed')
            cat(paste("Mutations are keept if",temp.col,"is greater or equal then:", trehsold, sep=' '), 
                file=logFile, append=TRUE, sep = "\n")    
            mut.all <- mut.all[mut.all[,temp.col] >= trehsold | is.na(mut.all[,temp.col]),  ]
         } else if (num.ind == 'h') {
            message('Everything higher then ',trehsold,' in column ', temp.col, ' will be removed')
            cat(paste("Mutations are keept if",temp.col,"is lower or equal then:", trehsold, sep=' '), 
                file=logFile, append=TRUE, sep = "\n")   
            mut.all <- mut.all[mut.all[,temp.col] <= trehsold | is.na(mut.all[,temp.col]),  ]
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

   
   annot.cols <- which(grepl('#',colnames(mut.all)) )
   samples.cols <- which(!grepl('#',colnames(mut.all)))
   remove.sam <- which(colSums(mut.all[,samples.cols], na.rm=T) == 0)
   if(length(remove.sam)>0){     
      cols.to.exclude <- which(colnames(mut.all) %in% names(remove.sam)) 
      mut.all <- mut.all[,-cols.to.exclude] 
      annot.cols <- which(grepl('#',colnames(mut.all)) )
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      #log file - snps numeric annot keept
      cat("\n------------------------------------------\n",  file=logFile, append=TRUE, sep = "\n")
      cat(paste("Samples with 0 mutations are excluded. These samples are:", paste0(names(remove.sam), collapse=', '), sep=' '), 
          file=logFile, append=TRUE, sep = "\n") 
   }
   rm(remove.sam)
   ################
   ################
   ################
   agg.pos <- which(agg.name == colnames(mut.all))

   cases.pos <-  which( colnames(mut.all) %in% cases.samp)
   contr.pos <-  which(colnames(mut.all) %in% controls.samp)
   mut.cases <- mut.all[, c(agg.pos, cases.pos)]
   mut.contr <- mut.all[, c(agg.pos, contr.pos) ]
   mut.cases[mut.cases==2] <- 1
   mut.contr[mut.contr==2] <- 1
   colnames(mut.cases)[1] <- 'AGG'
   colnames(mut.contr)[1] <- 'AGG'
   
   #mut.all.cases <- as.data.table(mut.all.cases)
   mut2.cases <- aggregate(.~AGG, data=mut.cases,FUN=sum, na.rm=TRUE, na.action=NULL)
   #mut.all.controls <- as.data.table(mut.all.controls)
   mut2.contr <- aggregate(.~AGG, data=mut.contr, FUN=sum, na.rm=TRUE, na.action=NULL)
   
   mut2.cases$sum <- apply(mut2.cases[,-1], 1, function (x) sum(x!=0, na.rm=T))
   mut2.contr$sum <- apply(mut2.contr[,-1], 1, function (x) sum(x!=0, na.rm=T))
   
   mut2.cases<- mut2.cases[,c('AGG', 'sum')]
   mut2.contr<- mut2.contr[,c('AGG', 'sum')]        
   
   mut2.cases$from <- 'cases'
   mut2.contr$from <- 'controls'
   
   mut2.cases$sum.cor <- (mut2.cases$sum/cases.num)*100
   mut2.contr$sum.cor <- (mut2.contr$sum/control.num)*100
   
   agg.all <- unique(c(mut.all[,agg.name]))
   #agg.all <- gsub('[();+:.>].*','', agg.all )
   #agg.all <- unique(agg.all)
   
   mut2  <- rbind(mut2.cases, mut2.contr)
   mut2$AGG <- factor(mut2$AGG, levels=agg.all)
   
   #rm(mut2.cases)
   #rm(mut2.contr)
   
   p2 <- ggplot(mut2, aes(sum.cor, fill=from)) + geom_density(alpha=0.5) + xlim(0,100) + theme_bw() 
   #save(p2, file='plot_density_cases_contr.RData')
   pdf('density_cases_contr.pdf', width=10, height=8)
   print(p2)
   dev.off()
   #log file
   cat(paste("density_cases_contr.pdf! - Density of mutations per gene, normalized to 100 samples - final"),
    file=logFile, append=TRUE, sep = "\n")

   mut3  <- merge(x=mut2.cases, y=mut2.contr, by='AGG',all=T)
   p3 <- ggplot(mut3, aes(x=sum.x, y=sum.y)) + geom_point(alpha=0.5) +  theme_bw() +
      geom_abline(intercept = 0, slope = control.num/cases.num, col=2) +
      xlab('num mut in genes in cases') + ylab('num mut in genes in controls')
   #save(p3, file='plot_counts_cases_contr.RData')
   
   pdf('counts_cases_contr.pdf' , width=10, height=8) 
   print(p3)
   dev.off()
   #log file
   cat(paste("counts_cases_contr.pdf! - Each dot is a gene, and x and y axis counts of mutations in cases and contorls - final"),
    file=logFile, append=TRUE, sep = "\n")
   
   mut3[is.na(mut3$sum.x),'sum.x'] <-0
   mut3[is.na(mut3$sum.cor.x),'sum.cor.x'] <-0
   mut3[is.na(mut3$sum.y),'sum.y'] <-0
   mut3[is.na(mut3$sum.cor.y),'sum.cor.y'] <-0
   
   cat(paste("Correlation betwen x and y axis (counts of mutations in cases and contorls ) is ", 
             cor(mut3$sum.x,mut3$sum.y)),
       file=logFile, append=TRUE, sep = "\n")
}
rm(answ)

cat("\nEverything done! Bye Bye ... \n", file=logFile, append=TRUE, sep = "\n")
cat("\nEverything done! Bye Bye ... \n")
