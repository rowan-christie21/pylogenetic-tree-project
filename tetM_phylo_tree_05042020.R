#5/2/2020
#Calculating CG and amino acid content of tetM sequences

#-------------------------------------------
# load libraries

library(ggplot2)
library(dplyr)
library(wesanderson)
library(RColorBrewer)
library(MuMIn)
library(car)
library(nlme)  
library(ncf)
library(lme4)
library(Matrix)
library(coefplot)
library(MASS)
library(genbankr)


#-------------------------------------------

tetm_strains <- readxl::read_xlsx("D:/Bioinformatics-Projects/phylogentic-tree-project-02292020/bacterial tetM datasets calculated 5_3_2020.xlsx", sheet=1)
tetm_strains <- tetm_strains[1:37,]

accs <- tetm_strains$`Accession number`

Accession <- c()
GC <- c()
seqs <- c()

for(i in accs) {
  id = GBAccession(i)
  #print(id)
  
  gb = readGenBank(id, partial = TRUE)
  gb_cds_tetM <- subset(gb@cds, gene=="tetM")
  
  start <- gb_cds_tetM@ranges@start
  end <- gb_cds_tetM@ranges@width+gb_cds_tetM@ranges@start-1
  seq <- gb@sequence
  cds_seq <- subseq(seq, start = start, end = end)
  
  Accession <- c(Accession, i)
  GC <- c(GC, ((alphabetFrequency(cds_seq)[, 2]/width(cds_seq))+(alphabetFrequency(cds_seq)[, 3]/width(cds_seq))))
  seqs <- c(seqs, cds_seq)
}

GCcont <- data.frame(Accession, GC)

seqlist<- DNAStringSetList(seqs)

tetm_tree <- msa(seqlist)
print(msa_res)
msaPrettyPrint(msa_res)
texi2pdf("msa_res.tex", clean=TRUE)

