######################################################################################
# CRG 
# Authors: Hana SUSAK, Georgia ESCARAMIS 
# date: 25/06/2018
#------------------------------------------------------------------------------------
# Correcting Gene names , to HUGO standards
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Calling the script: Rscript CorrectGeneNames.R --h
#------------------------------------------------------------------------------------
######################################################################################

rm(list=ls())
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", "--mut_file", type="character", help="input mutation file, QC output", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-c", "--col_num", type="integer", help="column nuber of gene names", metavar="number", required=TRUE)
parser$add_argument("-u", "--DB_username", type="character", help="username for databese", metavar="file", nargs=1)
parser$add_argument("-d", "--DB_password", type="character", help="password for database", metavar="file", nargs=1)
parser$add_argument("-s", "--syno", type="character", help="synonius mapping file", metavar="file", nargs=1)
parser$add_argument("-p", "--prev", type="character", help="previous mapping file", metavar="file", nargs=1)
parser$add_argument("-v", "--verbose", type="character", help="to print or not conversions (T or F)", metavar="boolean", nargs=1, default="F")


######################################################################################
# read input
######################################################################################
## get command line arguments
args <- parser$parse_args()

mut.file.name <- args$mut_file
# mut.file.name <- '/no_backup/GD/projects/CLL/germ_analysis/data/QC_output/mutations_filtered_samples.txt'

column.names <- as.integer(args$col_num)
# column.names <- 1

if(!is.null(args$DB_username) | !is.null(args$syno) ){
   internet.update <- FALSE
   from.db <- !is.null(args$DB_username)      
   from.file <- !is.null(args$syno) 
} else {
   internet.update <- TRUE
   print('you are updating previous and synonymous from link:  http://tinyurl.com/p8mh79k ')
   from.db <- FALSE
   from.file <- FALSE
}

print(paste('file name: ',mut.file.name))

print(paste('Column number: ',column.names))

ploting <- FALSE
printing <- as.logical(args$verbose) 
if(is.na(printing)){
  stop("-v (--verbose) input parameter can be only T or F")
}

if(internet.update){
   # choosen column from www.genenames.org:
   # Approved Symbol,   Locus Group,	Previous Symbols,	Synonyms
    url <-  'http://tinyurl.com/p8mh79k'
   df.mapping <- read.table(file=url, header=TRUE, fill=TRUE, sep='\t', na.strings='', stringsAsFactors=FALSE, quote="")
   #prev.df <- syn.prev[!is.na(syn.prev$Previous.Symbols),c('Approved.Symbol','Previous.Symbols'  )]
   #syno.df <- syn.prev[!is.na(syn.prev$Synonyms),c('Approved.Symbol','Synonyms'  )]
   apr.names <- df.mapping$Approved.Symbol
   df.mapping <- df.mapping[!(is.na(df.mapping$Previous.Symbols) & is.na(df.mapping$Synonyms)),]
   
   prev <- c()
   syno <- c()
   double.prev <- c()
   double.syno <- c()
   prev.approved <- c()
   syno.approved <- c()
   
   # mapping Previous nad Synonyms
   for(i in 1:nrow(df.mapping)){
      #previous names mapping
      if(! is.na(df.mapping[i,]$Previous.Symbols)){
         prev.names <- unlist(strsplit(df.mapping[i,]$Previous.Symbols, split=', '))
         for(prev.name in prev.names){
            if(! prev.name %in% apr.names){
               if (prev.name %in% names(prev)) {
                  double.prev <- c(double.prev, prev.name)
               }else {
                  prev[prev.name] <- df.mapping[i,1]
               } 
            } else {
               prev.approved <- c(prev.approved, prev.name)
            }
         }           
      }
      # synonim names
      if(! is.na(df.mapping[i,]$Synonyms)){
         syno.names <- unlist(strsplit(df.mapping[i,]$Synonyms, split=', '))
         for(syno.name in syno.names){
            if(! syno.name %in% apr.names){
               if (syno.name %in% names(syno)) {
                  double.syno <- c(double.syno, syno.name)
               } else {
                  syno[syno.name] <- df.mapping[i,1]
               }
            } else {
               syno.approved <- c(syno.approved, syno.name)
            }
         }
      }  
   }
   
   rm(i)
   
   double.prev <- unique(double.prev)
   double.syno <- unique(double.syno)
   prev.approved <- unique(prev.approved)
   syno.approved <- unique(syno.approved)
   
   prev <- prev[!names(prev) %in% double.prev]
   syno <- syno[!names(syno) %in% double.syno]
   
   #correct when one name have and previous and synonymus
   count <- 0
   names.in.both <- intersect(names(syno), names(prev))
   gene.mapping.problem <- c()
   
   for(gene in names.in.both){
      if (printing) { print('******************************************************') }
      if (printing) { print(paste('Gene is: ', gene )) }
      
      if (printing) { print('-----------------------------------')}
      prev.ap1.name <- prev[[gene]]    
      prev.loc.gr.1 <- df.mapping[df.mapping$Approved.Symbol == prev.ap1.name,]$Locus.Group
      if (printing) { print(paste('If previous name: ', prev.ap1.name)) }
      if (printing) { print(paste('Locus group: ', prev.loc.gr.1)) }
      
      if (printing) { print('-----------------------------------')}
      syno.ap2.name <- syno[[gene]]
      syno.loc.gr.2 <- df.mapping[df.mapping$Approved.Symbol == syno.ap2.name,]$Locus.Group
      if (printing) { print(paste('If synonimous name: ', syno.ap2.name)) }
      if (printing) { print(paste('Locus group: ', syno.loc.gr.2)) }
      
      if(prev.loc.gr.1==syno.loc.gr.2){
         prev <- prev[!names(prev)==gene]
         syno <- syno[!names(syno)==gene]
         gene.mapping.problem <- c(gene.mapping.problem, gene)
         count <- count + 1
      } else if (prev.loc.gr.1 == "protein-coding gene") {
         syno <- syno[!names(syno)==gene]
      } else if (syno.loc.gr.2 == "protein-coding gene"){
         prev <- prev[!names(prev)==gene]
      } else {
         prev <- prev[!names(prev)==gene]
         syno <- syno[!names(syno)==gene]
         gene.mapping.problem <- c(gene.mapping.problem, gene)
         count <- count + 1
      }
      
   }
 
   #write.table(as.data.frame(prev), file='map_previous_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(syno), file='map_synonymous_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(apr.names), file='aproved_names.txt'), sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(gene.mapping.problem), file='syno_and_prev_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(double.prev), file='double_prev_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(double.syno), file='double_syno_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(prev.approved), file='prev_approved_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)
   #write.table(as.data.frame(syno.approved), file='syno_approved_names.txt', sep='\t',  row.names=T, col.names=F,  quote=F)

   problem.genes <- unique(c(gene.mapping.problem, double.prev,double.syno , prev.approved, syno.approved))

   rm(count)
   rm(gene)
   rm(prev.ap1.name)
   rm(prev.loc.gr.1)
   rm(syno.ap2.name)
   rm(syno.loc.gr.2)
   rm(names.in.both)
   #rm(df.mapping)
   rm(url)
   rm(syno.names)
   rm(syno.name)
   rm(prev.names)
   rm(prev.name)
   rm(apr.names)
}

if(from.db ){
   suppressPackageStartupMessages(library(DBI,lib.loc='/software/so/el6.3/R-Modules'))
   suppressPackageStartupMessages(library(RMySQL, lib.loc='/software/so/el6.3/R-Modules'))

   print('You are checking pervious and synonymous from DB')
   print(paste('DB user name: ',DB_username))
   
   password.given <- DB_password
   
   query.mapping <- paste('SELECT gene_symbol, previous_symbol, Synonyms, LocusGroup FROM eDiVa_innoDB.Table_Gene_symbol_HGNC;')
   con <- dbConnect(MySQL(), user=user.given, password=password.given, dbname="eDiVa_innoDB", host="www.ediva.crg.eu")
   
   df.mapping  <- dbGetQuery(con, query.mapping) 
   
   dbDisconnect(con)
   
   apr.names <- df.mapping$gene_symbol
   
   #########
   prev <- c()
   syno <- c()
   
   
   for(i in 1:nrow(df.mapping)){
      #previous names mapping
      if(! is.na(df.mapping[i,]$previous_symbol)){
         prev.names <- unlist(strsplit(df.mapping[i,]$previous_symbol, split=', '))
         for(prev.name in prev.names){
            if(! prev.name %in% apr.names){
               prev[prev.name] <- df.mapping[i,1]
            }
         }
      }
      # synonim names
      if(! is.na(df.mapping[i,]$Synonyms)){
         syno.names <- unlist(strsplit(df.mapping[i,]$Synonyms, split=', '))
         for(syno.name in syno.names){
            if(! syno.name %in% apr.names){
               syno[syno.name] <- df.mapping[i,1]
            }
         }
      }  
   }
   
   
   #correct when one name have and previous and synonymus
   count <- 0
   names.in.both <- intersect(names(syno), names(prev))
   gene.mapping.problem <- c()
   
   for(gene in names.in.both){
      print('******************************************************')
      print(paste('Gene is: ', gene ))
      
      print('-----------------------------------')
      prev.ap1.name <- prev[[gene]]    
      prev.loc.gr.1 <- df.mapping[df.mapping$gene_symbol == prev.ap1.name,]$LocusGroup
      print(paste('If previous name: ', prev.ap1.name))
      print(paste('Locus group: ', prev.loc.gr.1))
      
      print('-----------------------------------')
      syno.ap2.name <- syno[[gene]]
      syno.loc.gr.2 <- df.mapping[df.mapping$gene_symbol == syno.ap2.name,]$LocusGroup
      print(paste('If synonimous name: ', syno.ap2.name))
      print(paste('Locus group: ', syno.loc.gr.2))
      
      if(prev.loc.gr.1==syno.loc.gr.2){
         prev <- prev[!names(prev)==gene]
         syno <- syno[!names(syno)==gene]
         gene.mapping.problem <- c(gene.mapping.problem, gene)
         count <- count + 1
      } else if (prev.loc.gr.1 == "protein-coding gene") {
         syno <- syno[!names(syno)==gene]
      } else if (syno.loc.gr.2 == "protein-coding gene"){
         prev <- prev[!names(prev)==gene]
      } else {
         prev <- prev[!names(prev)==gene]
         syno <- syno[!names(syno)==gene]
         gene.mapping.problem <- c(gene.mapping.problem, gene)
         count <- count + 1
      }
   }

}

#From  file
#first argument previous names file, second argument synonymous file name
if(from.file){
   file.name1 <- args$prev 
   print(paste('file name for synonimous names mapping: ',file.name1))
   
   file.name2 <- args$syno
   print(paste('file name for previous names mapping: ',file.name2))
   
   maps1 <- read.table(file=file.name1, header = F, sep='\t', stringsAsFactors=F,  strip.white = T, comment.char="")
   maps2 <- read.table(file=file.name2, header = F, sep='\t', stringsAsFactors=F,  strip.white = T,   quote = "", comment.char="")
   prev <- maps1[,2]
   names(prev) <-maps1[,1]
   syno <- maps2[,2]
   names(syno) <-maps2[,1]
   
   rm(maps1)
   rm(maps2)
}


correct.names <- function(x, map1=syno, map2=prev){
   gene.names.cor <- c()
   ex.gene <- c()
   for (gene in x){
      if (regexpr("[(]",gene)[1] > 0  ){
         if(printing){
            print(paste(gene, substr(gene, start=1, stop=(regexpr("[(]",gene)[1]-1))  , sep=' ----> '))
         }
         gene <- substr(gene, start=1, stop=(regexpr("[(]",gene)[1]-1))  
      } else if ( regexpr("[>]",gene)[1] > 0 ) {
         ex.gene <- c(ex.gene, gene) 
         if(printing){
            print(paste(gene,'wrong', sep=' ----> '))
         }
      } 
      if ( regexpr(";",gene)[1] > 0 ) {
         gene.temp <- gene
         gene1 <- substr(gene, start=1, stop=(regexpr("[;]",gene)[1]-1))
         gene2 <- substr(gene,  start=(regexpr("[;]",gene)[1]+1), stop=nchar(gene) )
         if (nchar(gene1) <= nchar(gene2)){
            gene.short <- gene1
            gene.long <- gene2
         } else {
            gene.short <- gene2
            gene.long <- gene1
         }    
         if (grepl(paste('-',gene.short, sep=""),gene.long) | grepl( paste(gene.short, '-', sep="" ),gene.long)  | grepl(gene.short ,gene.long )  ){
            gene <- gene.short
            if(printing){
               print(paste(gene.temp, gene  , sep=' ----> ')) 
            }
         } else {
            gene <- gene.temp
            if(printing){
               print(paste(gene.temp, gene  , sep=' --  TWO GENES --> '))  
            }
         }
      }
      if (gene %in% names(map1)){
         if(printing){
            print(paste(gene,map1[[gene]], sep=' -- synon --> '))
         }
         gene <- map1[[gene]]       
      }
      if (gene %in% names(map2)){
         if(printing){
            print(paste(gene,map2[[gene]], sep=' -- prev --> '))
         }
         gene <- map2[[gene]]
      }
      gene.names.cor <- c(gene.names.cor, gene)
   }
   return(list(new.names=gene.names.cor, exclude.names=ex.gene))
}


if(FALSE){
   correct.names <-  function(x, map1=syno, map2=prev){
      gene.names.cor <- c()
      ex.gene <- c()
      for (gene in x){
         if (regexpr("[(]",gene)[1] > 0  ){
            print(paste(gene, substr(gene, start=1, stop=(regexpr("[(]",gene)[1]-1))  , sep=' ---->'))
            gene <- substr(gene, start=1, stop=(regexpr("[(]",gene)[1]-1))  
         } else if ( regexpr("[>]",gene)[1] > 0 ) {
            ex.gene <- c(ex.gene, gene) 
            print(paste(gene,'wrong', sep=' ---->'))
         } 
         if ( regexpr(";",gene)[1] > 0 ) {
            gene.temp <- gene
            gene1 <- substr(gene, start=1, stop=(regexpr("[;]",gene)[1]-1))
            gene2 <- substr(gene,  start=(regexpr("[;]",gene)[1]+1), stop=nchar(gene) )
            if (nchar(gene1) <= nchar(gene2)){
               gene.short <- gene1
               gene.long <- gene2
            } else {
               gene.short <- gene2
               gene.long <- gene1
            }    
            if (grepl(paste('-',gene.short, sep=""),gene.long) | grepl( paste(gene.short, '-', sep="" ),gene.long)  | grepl(gene.short ,gene.long )  ){
               gene <- gene.short
               print(paste(gene.temp, gene  , sep=' ---->'))  
            } else {
               gene <- gene.temp
               print(paste(gene.temp, gene  , sep=' --  TWO GENES -->'))  
            }
         }
         if (gene %in% names(map1)){
            print(paste(gene,map1[[gene]], sep=' -- synon -->'))
            gene <- map1[[gene]]       
         }
         if (gene %in% names(map2)){
            print(paste(gene,map2[[gene]], sep=' -- prev -->'))
            gene <- map2[[gene]]
         }
         gene.names.cor <- c(gene.names.cor, gene)
      }
      return(list(new.names=gene.names.cor, exclude.names=ex.gene))
   }
   
}


if (substr(mut.file.name,nchar(mut.file.name)-3, nchar(mut.file.name)) == ".csv") {
   df <- read.csv(mut.file.name, sep=',', stringsAsFactors=F, quote="")
} else {
   df <- read.table(mut.file.name, sep='\t', stringsAsFactors=F, quote="", check.names=F, comment.char="", header=T)   
   xx <-  colnames(df)[column.names]
   colnames(df)[column.names] <- gsub('[;(+:.>-].*','', xx)
}

new.names2 <- correct.names(df[,column.names])
df[,column.names] <- new.names2$new.names
#df <- df[!df[,column.names] %in% new.names2$exclude.names, ]

if(grepl('.txt', mut.file.name)){
  new.file.name <- sub('.txt','_fix.txt', mut.file.name )
} else {
  new.file.name <- paste(mut.file.name,'fix', sep='_' )  
}



if (substr(mut.file.name,nchar(mut.file.name)-3, nchar(mut.file.name)) == ".csv"){
   write.table(df,file=new.file.name,  quote=F, sep=',', row.names=F,col.names = F )
} else {
   write.table(df,file=new.file.name,  quote=F, sep='\t', row.names=F, col.names = T )
}

print('Your file is created with new prefix!')
