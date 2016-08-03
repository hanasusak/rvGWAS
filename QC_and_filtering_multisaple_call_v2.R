######################################################################################
# CRG 
# Hana SUSAK
# date: 03/08/2015
#------------------------------------------------------------------------------------
# QC check for multisample call (PCA, filtering on sample info file, #mut per sample)
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Calling the script: Rscript QC_and_filtering_multisaple_call_v2.R --h
#------------------------------------------------------------------------------------
######################################################################################

rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(argparse))
#suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser()

parser$add_argument("-s", "--snp_file", type="character", help="input snp file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-i", "--indel_file", type="character", help="input indel file", metavar="file", nargs=1, required=FALSE)
parser$add_argument("-d", "--sample_desc_file", type="character", help="samples description file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-c", "--correct_sample_names", type="character", help="file with corrected names", metavar="file", nargs=1)

if (FALSE) {
   args <- commandArgs(trailingOnly = TRUE)
}

read.from <- 'stdin'
#read.from <- ''

######################################################################################
# read input
######################################################################################
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
snps.all <- read.table(file=snp.file.name, header=T, comment.char="", quote="", sep='\t', # nrows=100000,
                       stringsAsFactors=F, check.names=F)
cat("\nIn input file there are", nrow(snps.all), "unique snps  \n ....." )

refseq.pos <- which(grepl('\\(Refseq\\)',colnames(snps.all)))
new.names <- gsub('\\(Refseq\\)','',colnames(snps.all)[refseq.pos])
colnames(snps.all)[refseq.pos] <- new.names

annot.names.keep <- c('#Chr','#Position','#Reference','#Alteration','#Function','#Gene','#ExonicFunction',
                      '#EurEVSFrequency','#AfrEVSFrequency','#TotalEVSFrequency','#Eur1000GenomesFrequency',
                      '#Afr1000GenomesFrequency','#Asia1000GenomesFrequency','#Amr1000GenomesFrequency',
                      '#Total1000GenomesFrequency','#SegMentDup','#Cadd2')

if(any(!annot.names.keep %in% colnames(snps.all)) ){
   message(paste('You should have this column names in snps/indels file:', paste(annot.names.keep, collapse=', ')))
   stop('There are missing or wrongly named columns in your file')
}

annot.cols <- match(annot.names.keep, colnames(snps.all))
samples.cols <- which(!grepl('#',colnames(snps.all)))

snps.all <- snps.all[,c(annot.cols, samples.cols)]

if(!is.null(indel.file.name)){
   # InDels
   indels.all <- read.table(file=indel.file.name, header=T, comment.char="", quote="", sep='\t',  # nrows=100000, 
                            stringsAsFactors=F, check.names=F)
   cat("\nIn input file there are", nrow(indels.all), "unique indels \n ....." )
   
   refseq.pos <- which(grepl('\\(Refseq\\)',colnames(indels.all)))
   new.names <- gsub('\\(Refseq\\)','',colnames(indels.all)[refseq.pos])
   colnames(indels.all)[refseq.pos] <- new.names
     
   if(any(!annot.names.keep %in% colnames(indels.all)) ){
      message(paste('You should have this column names in snps/indels file:', paste(annot.names.keep, collapse=', ')))
      stop('There are missing or wrongly named columns in your file')
   }
   
   annot.cols <- match(annot.names.keep, colnames(indels.all))
   samples.cols <- which(!grepl('#',colnames(indels.all)))
   
   indels.all <- indels.all[,c(annot.cols, samples.cols)]
   
   #merge data sets
   indels.all$`#type` <- 'indel'
   snps.all$`#type` <- 'snp'
   
   if (ncol(indels.all)==ncol(snps.all) & all(colnames(snps.all) == colnames(indels.all))){
      mut.all <- rbind(snps.all, indels.all)
      
      annot.cols <- which(grepl('#',colnames(mut.all)))
      samples.cols <- which(!grepl('#',colnames(mut.all)))
      
      mut.all <- mut.all[,c(annot.cols, samples.cols)]
      
      rm(snps.all)
      rm(indels.all) 
   } else {
      stop('Column names in snps and indels file are not matching!')
   } 
} else {
   snps.all$`#type` <- 'snp'
   mut.all <- snps.all
   mut.all[mut.all$`#ExonicFunction` %in% c('frameshift deletion', 'frameshift insertion',
                                            'frameshift substitution','nonframeshift deletion',
                                            'nonframeshift insertion' , 'nonframeshift substitution') &
              !is.na(mut.all$`#ExonicFunction`), '#type'] <-  'indel'
   
   
   annot.cols <- which(grepl('#',colnames(mut.all)))
   samples.cols <- which(!grepl('#',colnames(mut.all)))   
   mut.all <- mut.all[,c(annot.cols, samples.cols)]
  
   rm(snps.all)   
}

rm(refseq.pos)
rm(new.names)

#mut.all <- mut.all[!duplicated(mut.all),]
#cat("\nThere are", nrow(mut.all), "unique snps/indels \n ....." )

mut.all <- mut.all[mut.all$`#Function` %in% c('exonic','exonic;splicing','splicing') | is.na(mut.all$`#Function`), ]
cat("\nAfter keeping only exonic and splicint there are", nrow(mut.all), "unique snps/indels \n ....." )

mut.all <- mut.all[mut.all$`#SegMentDup` < 0.9 | is.na(mut.all$`#SegMentDup`), ] 
cat("\nAfter filtering Segmental Duplication regions there are", nrow(mut.all), "unique snps/indels \n ....." )

annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

######################################################################################
# correcting some names
######################################################################################

if(exists('correct.names')){
   # correct some snps names
   snps.all.colnames <- colnames(mut.all)
   cor.index <- which(snps.all.colnames %in% rownames(correct.names))
   cor.name <- correct.names[ snps.all.colnames[snps.all.colnames %in% rownames(correct.names)] , 2]
   
   if (length(cor.index) > 0){
      cat("\nSample names in files: ", paste(colnames(mut.all)[cor.index], collapse=', '), 
          '\nare corrected to:',paste(cor.name, collapse=', ')," \n .....")
      colnames(mut.all)[cor.index] <- cor.name  
   }
   rm(cor.index)
   rm(cor.name)
   rm(snps.all.colnames)
   rm(correct.names)
}


######################################################################################
# filtering by info file
######################################################################################

samples.to.exclude <- c()

# collect kit info for samples
if(exists('samples.info.file')){
   samples.info <- read.table(file=samples.info.file, header=T, sep="\t", quote="", fill=T)
   message("Samples info/description file must contain samples IDs")
   cat("\nWhat is column with samples id (like in multisample snps/indels file) :", 
       paste('\n',1:ncol(samples.info),'.',colnames(samples.info), sep=''),
       "\nChoose number: ")
   answ.col <- scan(read.from,what=integer(),n=1,quiet=TRUE)
   answ.col <- answ.col[answ.col<= ncol(samples.info)]       
   answ.name <- colnames(samples.info)[answ.col]
   cat("\nYou chose that samples id column is", answ.name, " ... ") 
   answ <- 'y'
   #colnames(samples.info)[answ.col] <- 'samples'
   if( any(!colnames(mut.all)[samples.cols] %in% samples.info[,answ.col])){
      stop(' In samples info file there must be all samples from mutation file')
   } else {
      samples.info <- samples.info[samples.info[,answ.col] %in%  colnames(mut.all)[samples.cols],]      
   }
   
   cat("\n\nYou can filter by this colums :", 
       paste('\n',1:ncol(samples.info),'.',colnames(samples.info), sep=''),
       "\nYou want to filer (y or n): ")
   answ <- scan(read.from,what=character(),n=1,quiet=TRUE)
   if(answ=='y'){
      cat("\nEnter column numbers space separated: ")
      x <- scan(read.from,what=integer(),nmax=ncol(samples.info),quiet=TRUE, nlines=1)
      x <- unique(x)
      x <- x[x<= ncol(samples.info)]       
      if (length(x) <= 0){
         cat("\nYou or chose not to filter or wrong columns \n.....")
      } else {
         cat("\nColumn(s) you chose is:",paste(colnames(samples.info)[x], collapse=', ' ), "\n.....")
         for (i in x){
            options <- as.character(unique(samples.info[,i]))
            options[is.na(options)] <- 'NA_Missing_info'
            cat("\n\nIn column",colnames(samples.info)[i],"there are folowing fields:")
            counts.tb <- table(samples.info[,i], useNA='always')
            names(counts.tb)[is.na(names(counts.tb))] <- 'NA_Missing_info'
            counts.tb <- counts.tb[options]
            cat('\n',paste(1:length(options),'.',(options),' (',counts.tb,' samples)','\n', sep=''))
            cat("\nEnter  numbers which are in front of options you don't want to keep (space separated): ")
            exclude.options <- scan(read.from,what=integer(),nmax=length(options),quiet=TRUE, nlines=1)
            exclude.options <- unique(exclude.options)
            exclude.options <- exclude.options[exclude.options <= length(options)]
            if (length(exclude.options) > 0) {
               cat("\nSamples whos", colnames(samples.info)[i],
                   "is in:", paste(options[exclude.options], collapse=', '),
                   "- will be excluded")  
               if( 'NA_Missing_info' %in% options[exclude.options]){
                  samples.to.exclude.temp <-  as.character(samples.info[is.na(as.character(samples.info[,i])) ,answ.name ])
                  samples.to.exclude <- unique(c(samples.to.exclude, samples.to.exclude.temp)) 
               }
               samples.to.exclude.temp <-  as.character(samples.info[as.character(samples.info[,i]) %in% options[exclude.options],answ.name ])
               samples.to.exclude <- unique(c(samples.to.exclude, samples.to.exclude.temp)) 
            }   else {
               cat("\n Nothing from column",colnames(samples.info)[i],"will be excluded ")
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
}


cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
cat("\nTotal number of samples to exclude is",length(cols.to.exclude)," \n.....")

if (length(cols.to.exclude) > 0){
   mut.all <- mut.all[,-cols.to.exclude]
   samples.info2 <- samples.info[!samples.info[,answ.col] %in% samples.to.exclude,]
   samples.info2 <- droplevels(samples.info2)
} else {
   samples.info2 <- samples.info
}
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

rm(samples.info)
rm(snp.file.name)
rm(indel.file.name)
rm(samples.info.file)
rm(sample.correct.file.name)

######################################################################################
# Filter outliers by number of mutations
######################################################################################

mut.all$`#snpID` <- apply(mut.all[,c("#Chr","#Position","#Reference","#Alteration")], 1, function(x) paste(trimws(x),collapse=';'))
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

# reorder, first annotation, then samples columns
mut.all <- mut.all[,c(annot.cols, samples.cols )]
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

cc <- which(colnames(mut.all)%in% c( "#ExonicFunction", "#type", "#snpID"))

for(i in 1:ceiling(nrow(mut.all)/100000)){
   max.row <- min((i*100000),nrow(mut.all)) 
   #print(max.row)
   melt.mut.all.temp <- melt(mut.all[((i-1)*100000+1):max.row,c(cc,samples.cols)], 
                             id.vars=c( "#ExonicFunction", "#type", "#snpID"))
   melt.mut.all.temp <- melt.mut.all.temp[melt.mut.all.temp$value != 0,]
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

all.samples <- as.character(colnames(mut.all[,samples.cols]))
melt.mut.all$variable <- factor(melt.mut.all$variable, levels=all.samples)

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

X11(width=14, height=10)
p1 <- ggplot(samples.info2, aes(x=Num.Mutations)) + geom_histogram() + theme_bw()
p1

cat("\nEnter minimum number of mutations per sample: ")
min.mut <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter maximum number of mutations per sample: ")
max.mut <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)

dev.off()
#plot hist of per saple number of mutations
pdf(file='hist.number.mutations_old1.pdf', width=15,height=12,useDingbats=F)
p1
dev.off()
rm(p1)

samples.to.exclude.temp <- c()
if (length(min.mut) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(samples.info2[samples.info2$Num.Mutations < min.mut,answ.name]))
}
if (length(max.mut) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(samples.info2[samples.info2$Num.Mutations > max.mut,answ.name]))
}
samples.to.exclude.temp <- unique(samples.to.exclude.temp)

cat("Samples you decided to exclude because of number of mutations are:",
    paste(samples.to.exclude.temp, collapse=", "),"\n.....")

samples.to.exclude <- c(samples.to.exclude, samples.to.exclude.temp)
rm(samples.to.exclude.temp)
rm(max.mut)
rm(min.mut)

cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
cat("\nTotal number of samples to exclude because of number of mutations is",length(cols.to.exclude),"\n.....")

if (length(cols.to.exclude) > 0){
   mut.all <- mut.all[,-cols.to.exclude]   
   melt.mut.all <- melt.mut.all[!melt.mut.all$variable %in% samples.to.exclude,]
   melt.mut.all <- droplevels(melt.mut.all)
   samples.info2 <- samples.info2[!samples.info2[,answ.col] %in% samples.to.exclude,]
   samples.info2 <- droplevels(samples.info2)
   #all.samples <- all.samples[!all.samples %in% samples.to.exclude]
}

#recaluclate annotation and sample columns
annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))


pdf(file='hist.number.mutations_old2.pdf', width=15,height=12,useDingbats=F)
p1 <- ggplot(samples.info2, aes(x=Num.Mutations)) + geom_histogram() + theme_bw()
p1
dev.off()
rm(p1)

#################################################################################
# to plot in the end samples number of mutations distribution (colored by type)
#################################################################################

ordered.samples <- as.character(names(sort(table(melt.mut.all$variable), decreasing=T)))

melt.mut.all$type2 <- 'SNP_nonsilent'
melt.mut.all[melt.mut.all$`#ExonicFunction` == 'synonymous SNV' & !is.na(melt.mut.all$`#ExonicFunction`), 'type2'] <- 'SNP_silent'   
 
melt.mut.all[melt.mut.all$`#type` == 'indel', 'type2'] <- 'InDel'

melt.mut.all <- melt.mut.all[,c('variable', 'type2')]
colnames(melt.mut.all) <- c('Sample', 'Type')
temp.df <- table(melt.mut.all[,c('Sample','Type')])
melt.fin.df <- melt(temp.df)

melt.fin.df$Sample <- factor(melt.fin.df$Sample, levels=ordered.samples)

rm(temp.df)
rm(melt.mut.all)

samples.info2[,answ.col] <- factor(samples.info2[,answ.col], levels=ordered.samples)


######################################################################################
# PCA QC and filtering
######################################################################################
cat("\n\nYou  want to do PCA only on silent mutations? \nChoose y if yes, or n to do it on all mutations(y or n): ")
answ.pca <- scan(read.from,what=character(),n=1,quiet=TRUE)
if(answ.pca=='y'){ 
   silent.mut <- mut.all[ mut.all$`#ExonicFunction`=='synonymous SNV' & 
                             !is.na(mut.all$`#ExonicFunction`), samples.cols]
   
   if(nrow(silent.mut) > 500000){
      silent.mut <- silent.mut[sample(nrow(silent.mut),500000),]
      cat("\nPCA is caluclated on random subeset of 200K silent mutations ")
   } else {
      cat("\nPCA is caluclated on silent mutaions! ")  
   }
} else {
   #not realy silent, just keep name for siplicity
   silent.mut <- mut.all[, samples.cols]
   if(nrow(silent.mut) > 500000){
      silent.mut <- silent.mut[sample(nrow(silent.mut),500000),]
      cat("\nPCA is caluclated on random subeset of 500K all mutations ")
   } else {
      cat("\nPCA is caluclated on all mutaions! ")  
   }
}   


mat.pat <- t(as.matrix(silent.mut))

rm(silent.mut)
# exclud patients without any mutation
#mat.pat <- mat.pat[rowSums(mat.pat) != 0, ]

# exclude snps with no patient mutated
mat.pat <- mat.pat[ ,colSums(mat.pat) != 0]
var.nip <- tune.pca(X=mat.pat, center=F )
dev.off()

pc.max <- min(dim(mat.pat))
X11(width=14, height=10)
par(mfrow=c(2,1))
plot(var.nip$cum.var, type='l', ylim=c(0,1), col=2, xlab="PC components", ylab="Proportion of Explaine Variance")
plot(var.nip$cum.var, type='o', pch=19, xlim=c(0,50), 
     ylim=c(min(var.nip$cum.var[1:min(50,pc.max)])*0.98, max(var.nip$cum.var[1:min(50,pc.max)])*1.02),
     xlab="PC components (first 50)", ylab="Proportion of Explaine Variance")

# ask user how many components tey want
cat("\nEnter  numbers of principal components you want to be calculated (n <= 50): \n (Warrning: only first two PC comp. will be plotted for filtering) ")
n.cp <- scan(read.from,what=integer(),nmax=1,quiet=TRUE, nlines=1)
n.cp <- min(c(n.cp, 50))
dev.off()

# first 2 copmponents of PCA
mat.pat.nip <- nipals(X=mat.pat, ncomp=2, reconst=T)

df <- data.frame( samples.pca = rownames(mat.pat),
                  PC1 = mat.pat.nip$t[,1],
                  PC2 = mat.pat.nip$t[,2])


pca.snps.all <- merge(x=df, y=samples.info2, by.x='samples.pca', by.y=answ.name,  all.x=T)
rm(df)
rm(mat.pat)

# choose color and shape columns
cat("\n\nYou can color/shape in PCA plot by:", 
    paste('\n',1:(ncol(pca.snps.all)-3),'.',colnames(pca.snps.all)[-(1:3)], sep=''))
cat("\nEnter  numbers space separated for color and shape (in this order, or one/none if you don't want color/shape): ")
x <- scan(read.from,what=integer(),nmax=ncol(2),quiet=TRUE, nlines=1)
if (length(x) > 0){
   color.col <- colnames(pca.snps.all[,-c(1:3)])[x[1]]
} 

if (length(x) > 1){
   shape.col <- colnames(pca.snps.all[,-c(1:3)])[x[2]]
} 

if (exists('shape.col')){
   num.shapes <- length(unique(pca.snps.all[,x[2]+3]))
   if (num.shapes > 25 ){
      vv <- c(sample(1:25),sample(1:25,size=(num.shapes-25), replace=F))   
   } else {
      vv <- c(sample(1:num.shapes))
   }
   rm(num.shapes)
}
rm(x)

######

cat("\n\nIn the plot choose cordinates in which you want to keep the data. \n
    (is there data points which look like outliers? Choose thresholds for x-min, x-max, y-min and y-max) \n .....")


X11(width=15, height=10)
if (exists('color.col') & exists('shape.col')){
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', colour=color.col, shape=shape.col )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
   p2
} else if (exists('color.col')) {
   p2 <- ggplot(data=pca.snps.all2, aes_string(x='PC1', y='PC2', colour=color.col )) +
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

#p3 <-ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC3', colour=color.col, shape=shape.col )) +
#   geom_point(size=4, alpha=0.8) +
#   scale_shape_manual(values=vv) + theme_bw()
#p33 <-ggplot(data=pca.snps.all2, aes_string(x='PC1', y='PC3', colour=color.col, shape=shape.col )) +
#   geom_point(size=4, alpha=0.8) +
#   scale_shape_manual(values=vv2) + theme_bw()



cat("\nEnter x-min  value (just press enter if you don't want any): ")
xmin <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter x-max value(just press enter if you don't want any): ")
xmax <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter y-min value(just press enter if you don't want any): ")
ymin <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)
cat("\nEnter y-max value(just press enter if you don't want any): ")
ymax <- scan(read.from,what=numeric(),nmax=ncol(1),quiet=TRUE, nlines=1)


dev.off()

pdf(file='projects.pca.allmut_old1.pdf', width=20,height=20,useDingbats=F)
p2
dev.off()
rm(p2)

#pdf(file='projects.pca.cancers_old3.pdf', width=20,height=20,useDingbats=F)
#p33
#dev.off()

samples.to.exclude.temp <- c()
if (length(xmin) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC1 < xmin,'samples.pca']))
}
if (length(xmax) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC1 > xmax,'samples.pca']))
}
if (length(ymin) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC2 < ymin,'samples.pca']))
}
if (length(ymax) >0 ){
   samples.to.exclude.temp <- c(samples.to.exclude.temp,as.character(pca.snps.all[pca.snps.all$PC2 > ymax,'samples.pca']))
}
samples.to.exclude.temp <- unique(samples.to.exclude.temp)

cat("Samples (",length(samples.to.exclude.temp)," of them ) you decided to exclude because of PCA plot are:",paste(samples.to.exclude.temp, collapse=", "), "\n .....")
samples.to.exclude <- c(samples.to.exclude, samples.to.exclude.temp)
rm(samples.to.exclude.temp)


cols.to.exclude <- which(colnames(mut.all) %in% samples.to.exclude)
cat("\nTotal number of samples to exclude because of PCA is",length(cols.to.exclude), "\n .....")


if (length(cols.to.exclude) > 0){
   mut.all <- mut.all[,-cols.to.exclude]
   
   melt.fin.df <- melt.fin.df[!melt.fin.df$Sample %in% samples.to.exclude,]
   melt.fin.df$Sample <- droplevels(melt.fin.df$Sample)
   
   samples.info2 <- samples.info2[!samples.info2[,answ.col] %in% samples.to.exclude,]
   samples.info2 <- droplevels(samples.info2)
   
   #ordered.samples <- ordered.samples[!ordered.samples %in% samples.to.exclude ]
   #samples.info2[,answ.col] <- factor(samples.info2[,answ.col], levels=ordered.samples)

}
rm(pca.snps.all)

annot.cols <- which(grepl('#',colnames(mut.all)) )
samples.cols <- which(!grepl('#',colnames(mut.all)))

######################################################################################
# Repolot everything after final filtration
######################################################################################
pdf(file='sample.mutation.distribution.pdf', width=25,height=15,useDingbats=F)
ggplot(data=melt.fin.df, aes(x=Sample, y=value, fill=Type)) + geom_bar(stat="identity")   + 
   theme_bw() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
   scale_y_continuous(expand=c(0,0)) + xlab('Samples') + ylab('Number of mutations')
dev.off()  
rm(melt.fin.df)


#plot hist of per saple number of mutations
pdf(file='hist.number.mutations.pdf', width=12,height=12,useDingbats=F)
ggplot(samples.info2, aes(x=Num.Mutations)) + geom_histogram() + theme_bw()
dev.off()


##############
# repeat PCA
##############
if(answ.pca=='y'){ 
   silent.mut <- mut.all[ mut.all$`#ExonicFunction`=='synonymous SNV' & 
                             !is.na(mut.all$`#ExonicFunction`), samples.cols]
   
   if(nrow(silent.mut) > 500000){
      silent.mut <- silent.mut[sample(nrow(silent.mut),500000),]
      cat("\nPCA is recaluclated on random subeset of 200K silent mutations ")
   } else {
      cat("\nPCA is recaluclated on silent mutaions! ")  
   }
} else {
   #not realy silent, just keep name for siplicity
   silent.mut <- mut.all[, samples.cols]
   if(nrow(silent.mut) > 500000){
      silent.mut <- silent.mut[sample(nrow(silent.mut),500000),]
      cat("\nPCA is recaluclated on random subeset of 500K all mutations ")
   } else {
      cat("\nPCA is recaluclated on all mutaions! ")  
   }
}  


mat.pat <- t(as.matrix(silent.mut))

rm(silent.mut)

# exclude snps with no patient mutated
mat.pat <- mat.pat[ ,colSums(mat.pat) != 0]
var.nip <- tune.pca(X=mat.pat, center=F)
dev.off()

pc.max <- min(dim(mat.pat))

pdf(file='PCA_variance_explained.pdf', width=12,height=12,useDingbats=F)
par(mfrow=c(2,1))
plot(var.nip$cum.var, type='l', ylim=c(0,1), col=2, xlab="PC components", ylab="Proportion of Explaine Variance")
#abline(v=n.cp, col="blue")
plot(var.nip$cum.var, type='o', pch=19, xlim=c(0,50), 
     ylim=c(min(var.nip$cum.var[1:min(50,pc.max)])*0.98, max(var.nip$cum.var[1:min(50,pc.max)])*1.02),
     xlab="PC components (first 50)", ylab="Proportion of Explaine Variance")
abline(v=n.cp, col="blue")
dev.off()

# first n copmponents of PCA
n.cp2 <- max(2,n.cp)
mat.pat.nip <- nipals(X=mat.pat, ncomp=n.cp2, reconst=T)

df <- cbind( data.frame(samples.pca= rownames(mat.pat)),mat.pat.nip$t)
colnames(df) <- c('samples.pca', paste('PC',1:ncol(mat.pat.nip$t), sep=""))

pca.snps.all <- merge(x=df, y=samples.info2,  by.x='samples.pca', by.y=answ.name, all.x=T)


#plot pca
if (exists('color.col') & exists('shape.col')){
   p2 <- ggplot(data=pca.snps.all, aes_string(x='PC1', y='PC2', colour=color.col, shape=shape.col )) +
      geom_point(size=4, alpha=0.8) +
      scale_shape_manual(values=vv) + theme_bw()
} else if (exists('color.col')) {
   p2 <- ggplot(data=pca.snps.all2, aes_string(x='PC1', y='PC2', colour=color.col )) +
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
pdf(file='projects.pca.allmut.pdf', width=18,height=15,useDingbats=F)
p2
dev.off()
rm(p2)

rm(mat.pat)
rm(mat.pat.nip)
rm(pc.max)
rm(xmax)
rm(xmin)
rm(ymax)
rm(ymin)

cat("\nThere is left ", nrow(mut.all)," mutations and indels together ")
cat("\nTable dimensions are: ", dim(mut.all), "(",ncol(mut.all)-length(annot.cols)," samples ) \n.....")

write.table(mut.all, file='mutations_filtered_samples.txt', quote=F, sep='\t', row.names=F)

# correct it only one PCA is chosen
if(n.cp2 != n.cp) {
   df <- df[,1:(n.cp+1)] 
   pca.snps.all <- merge(x=df, y=samples.info2,  by.x='samples.pca', by.y=answ.name, all.x=T)
}
rm(df)

write.table(pca.snps.all, file='samples_info_pca.txt', quote=F, sep='\t', row.names=F)


######################################################################################
# distribution cases vs controls
######################################################################################
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
   cat("\n\nIn column",colnames(samples.info2)[x],"there are folowing fields:")
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
          "- will be considered as samples")              
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
   samples.info2[ samples.info2[,answ.col] %in% cases ,'status'] <- 'cases'
   
   samples.info2[,answ.name] <- factor(samples.info2[,answ.name], 
                             levels=unique(as.character(samples.info2[order(-samples.info2$Num.Mutations) ,answ.name])))
   
   p1 <- ggplot(data=samples.info2, aes_string(x=answ.name, y='Num.Mutations', fill='status')) + geom_bar(stat="identity")   + 
      theme_bw() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
      scale_y_continuous(expand=c(0,0)) + xlab('Samples') + ylab('Number of mutations')
   #save(p1, file='plot_sample.mutation.hist_caseVScontrol.RData')
   
   pdf(file='sample.mutation.hist_caseVScontrol.pdf', width=25,height=15,useDingbats=F)
   print(p1)
   dev.off() 
 
   
}



if(answ=='y' & TRUE ){
   #density of mutations per gene
   cases.num <- sum(samples.info2$status == 'cases')
   #cases.samp <- as.character(samples.info2[samples.info2$status == 'cases','samples'])
   cases.samp <- cases
   control.num <-  sum(samples.info2$status == 'controls')
   controls.samp <- as.character(samples.info2[samples.info2$status == 'controls',answ.name])
   
   gene.pos <- which('#Gene' == colnames(mut.all))
   #ef.pos <- which('#ExonicFunction' == colnames(mut.all))
   cadd.pos <- which('#Cadd2' == colnames(mut.all))
   
   cases.pos <-  which( colnames(mut.all) %in% cases.samp)
   contr.pos <-  which(colnames(mut.all) %in% controls.samp)
   
   cat("\nDo you want now to check distribution per gene  on silent (s), nonsilent (n) or all (a)  - you can put space saprated several options: ")
   answ2 <- (scan(read.from,what=character(),nmax=3,nlines=1,quiet=TRUE))
   answ2 <- answ2[answ2 %in% c('a','s','n')]
   if (length(answ2)==0)  {
      stop("you didn't choose corecct options")
   }

   
   for (op in answ2){
      print(op)
      if ('s' == op  ){
         mut.cases <- mut.all[mut.all$`#ExonicFunction`  == 'synonymous SNV' & 
                                 !is.na(mut.all$`#ExonicFunction`), 
                              c(gene.pos, cadd.pos, cases.pos)]
         mut.contr <- mut.all[mut.all$`#ExonicFunction`  == 'synonymous SNV' &
                                 !is.na(mut.all$`#ExonicFunction`), 
                              c(gene.pos, cadd.pos, contr.pos) ]
         txt <- 'silent'
      } else if ( 'n' == op  ){
         mut.cases <- mut.all[mut.all$`#ExonicFunction`  != 'synonymous SNV' | 
                                 is.na(mut.all$`#ExonicFunction`), 
                              c(gene.pos, cadd.pos, cases.pos)]
         mut.contr <- mut.all[mut.all$`#ExonicFunction`  != 'synonymous SNV' |
                                 is.na(mut.all$`#ExonicFunction`), 
                              c(gene.pos, cadd.pos, contr.pos) ]
         txt <- 'nonsilent'
      } else  if ('a' == op ){
         mut.cases <- mut.all[, c(gene.pos,cadd.pos, cases.pos)]
         mut.contr <- mut.all[, c(gene.pos,cadd.pos, contr.pos) ]
         txt <- 'all'
      } 
   
      cat("\n--------------------------------------- ")
      cat(paste("\nCadd score threshold (0 will not filter) for ",txt," mutations: ", sep=""))
      cadd2 <- as.numeric(scan(read.from,what=character(),n=1,quiet=TRUE))
      if ( cadd2 > 0){
         mut.cases <- mut.cases[mut.cases$`#Cadd2` > cadd2 | is.na(mut.cases$`#Cadd2`), -2]
         mut.contr <- mut.contr[mut.contr$`#Cadd2` > cadd2 | is.na(mut.contr$`#Cadd2`), -2]
      } else {
         mut.cases <- mut.cases[, -2]
         mut.contr <- mut.contr[, -2]
      }
      
      
      colnames(mut.cases)[1] <- 'Gene'
      colnames(mut.contr)[1] <- 'Gene'
      mut.cases$Gene <- gsub('[();+:.>].*','', as.character(mut.cases$Gene) )
      mut.contr$Gene <- gsub('[();+:.>].*','', as.character(mut.contr$Gene) )
      # gsub('[+:(.;>-]+[A-z0-9()]+','', as.character(mut.contr$Gene)  )
      
      mut.cases[mut.cases==2] <- 1
      mut.contr[mut.contr==2] <- 1
   
      
      #mut.all.cases <- as.data.table(mut.all.cases)
      mut2.cases <- aggregate(.~Gene, data=mut.cases, sum)
      #mut.all.controls <- as.data.table(mut.all.controls)
      mut2.contr <- aggregate(.~Gene, data=mut.contr, sum) 
      
      mut2.cases$sum <- apply(mut2.cases[,-1], 1, function (x) sum(x!=0))
      mut2.contr$sum <- apply(mut2.contr[,-1], 1, function (x) sum(x!=0))
      
      mut2.cases<- mut2.cases[,c('Gene', 'sum')]
      mut2.contr<- mut2.contr[,c('Gene', 'sum')]
   
      
      mut2.cases$from <- 'cases'
      mut2.contr$from <- 'controls'
      
      mut2.cases$sum.cor <- (mut2.cases$sum/cases.num)*100
      mut2.contr$sum.cor <- (mut2.contr$sum/control.num)*100
      
      genes.all <- unique(c(mut.all$`#Gene`))
      genes.all <- gsub('[();+:.>].*','', genes.all )
      genes.all <- unique(genes.all)
      
      mut2  <- rbind(mut2.cases, mut2.contr)
      mut2$Gene <- factor(mut2$Gene, levels=genes.all)
      
      #rm(mut2.cases)
      #rm(mut2.contr)
   
      
      p2 <- ggplot(mut2, aes(sum.cor, fill=from)) + geom_density(alpha=0.5) + xlim(0,100) + theme_bw() 
      #save(p2, file='plot_density_cases_contr.RData')
      pdf(paste(txt,'density_cases_contr.pdf',sep="_"), width=10, height=8)
      print(p2)
      dev.off()
      
      mut3  <- merge(x=mut2.cases, y=mut2.contr, by='Gene',all=T)
      p3 <- ggplot(mut3, aes(x=sum.x, y=sum.y)) + geom_point(alpha=0.5) +  theme_bw() +
         geom_abline(intercept = 0, slope = control.num/cases.num, col=2) +
         xlab('num mut in genes in cases') + ylab('num mut in genes in controls')
      #save(p3, file='plot_counts_cases_contr.RData')
      
      pdf(paste(txt,'counts_cases_contr.pdf',sep="_") , width=10, height=8) 
      print(p3)
      dev.off()
   }
     
   
  
  rm(op)
  
  
   
}
rm(answ)


cat("\nEverything done! Bye Bye ... \n")
