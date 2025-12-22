####Package####

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(maftools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)

####WorkingDirection####
raw_WD <- '/public-SSD/home/sunpc/scseq/backend/R/TCGA'
setwd(raw_WD)

list.files()
if(!dir.exists('code')){
  dir.create('code')
  }

#####Data####

list.files('data',full.names = T)

TCGAlist <- list.files('data')
TCGAdir <- list.files('data',full.names = T)

TCGAinfo <- data.frame(Cancer=TCGAlist,
                       mRNA=rep(NA,length(TCGAlist)),
                       Methylation=rep(NA,length(TCGAlist)),
                       Protein=rep(NA,length(TCGAlist)),
                       Surv=rep(NA,length(TCGAlist)),
                       Clinical=rep(NA,length(TCGAlist)),
                       mRNAcount=rep(NA,length(TCGAlist)),
                       mRNAtpm=rep(NA,length(TCGAlist))
                       )

for (i in 1:length(TCGAlist)) {
  
  print(paste0('Accessing ',TCGAlist[i]))
  
  dir <- list.files(TCGAdir[i])
  
  #mRNA
  if(!identical(
    dir[grepl("HiSeqV2", dir) & !grepl(".gz", dir)],
    character(0)
    )
  ){
    TCGAinfo$mRNA[i] <- paste0(TCGAdir[i],'/',dir[grepl("HiSeqV2", dir) & !grepl(".gz", dir)])
    print(paste0(' mRNA: ',
                 paste0(TCGAdir[i],'/',dir[grepl("HiSeqV2", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$mRNA[i] <- 'none'
    print(" mRNA: No matching files found!")
  }
  
  #Methylation
  if(!identical(
    dir[grepl("HumanMethylation450", dir) & !grepl(".gz", dir)],
    character(0)
    )
  ){
    TCGAinfo$Methylation[i] <- paste0(TCGAdir[i],'/',dir[grepl("HumanMethylation450", dir) & !grepl(".gz", dir)])
    print(paste0(' Methylation: ',
                 paste0(TCGAdir[i],'/',dir[grepl("HumanMethylation450", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$Methylation[i] <- 'none'
    print(" Methylation: No matching files found!")
  }
  
  #Protein
  if(!identical(
    dir[grepl("RPPA", dir) & !grepl(".gz", dir)],
    character(0)
  )
  ){
    TCGAinfo$Protein[i] <- paste0(TCGAdir[i],'/',dir[grepl("RPPA", dir) & !grepl(".gz", dir)])
    print(paste0(' Protein: ',
                 paste0(TCGAdir[i],'/',dir[grepl("RPPA", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$Protein[i] <- 'none'
    print(" Protein: No matching files found!")
  }
  
  #Surv
  if(!identical(
    dir[grepl("survival", dir) & !grepl(".gz", dir)],
    character(0)
  )
  ){
    TCGAinfo$Surv[i] <- paste0(TCGAdir[i],'/',dir[grepl("survival", dir) & !grepl(".gz", dir)])
    print(paste0(' Surv: ',
                 paste0(TCGAdir[i],'/',dir[grepl("survival", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$Surv[i] <- 'none'
    print(" surv: No matching files found!")
  }
  
  #Clinical
  if(!identical(
    dir[grepl("clinicalMatrix", dir) & !grepl(".gz", dir)],
    character(0)
  )
  ){
    TCGAinfo$Clinical[i] <- paste0(TCGAdir[i],'/',dir[grepl("clinicalMatrix", dir) & !grepl(".gz", dir)])
    print(paste0(' Clinical: ',
                 paste0(TCGAdir[i],'/',dir[grepl("clinicalMatrix", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$Clinical[i] <- 'none'
    print(" Clinical: No matching files found!")
  }
  
  #count
  if(!identical(
    dir[grepl("star_counts", dir) & !grepl(".gz", dir)],
    character(0)
  )
  ){
    TCGAinfo$mRNAcount[i] <- paste0(TCGAdir[i],'/',dir[grepl("star_counts", dir) & !grepl(".gz", dir)])
    print(paste0(' mRNA_counts: ',
                 paste0(TCGAdir[i],'/',dir[grepl("star_counts", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$mRNAcount[i] <- 'none'
    print(" mRNA_counts: No matching files found!")
  }
  
  #TPM
  if(!identical(
    dir[grepl("star_tpm", dir) & !grepl(".gz", dir)],
    character(0)
  )
  ){
    TCGAinfo$mRNAtpm[i] <- paste0(TCGAdir[i],'/',dir[grepl("star_tpm", dir) & !grepl(".gz", dir)])
    print(paste0(' mRNA_TPM: ',
                 paste0(TCGAdir[i],'/',dir[grepl("star_tpm", dir) & !grepl(".gz", dir)])
    ))
  }else{
    TCGAinfo$mRNAtpm[i] <- 'none'
    print(" mRNA_TPM: No matching files found!")
  }
  
  print(paste0("Done: ",TCGAlist[i]))
}

####save####
if(!dir.exists('res')){
  dir.create('res')
}

write.csv(TCGAinfo,file = 'res/TCGAdir.csv',row.names = T)



