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

#read in bacterial data file
tetm_strains <- readxl::read_xlsx("D:/Bioinformatics-Projects/phylogentic-tree-project-02292020/bacterial tetM datasets 5_3_2020.xlsx", sheet=1)
#limit to records with MIC data
tetm_strains <- tetm_strains[1:37,]

#make list of accession numbers
accs <- tetm_strains$`Accession number`

#set up vectors
Accession <- c()
GC <- c()
seqs <- c()
aa_freq <- c()

A <- c()
C <- c()
D <- c()
E <- c()
F <- c()
G <- c()
H <- c()
I <- c()
K <- c()
L <- c()
M <- c()
N <- c()
P <- c()
Q <- c()
R <- c()
S <- c()
T <- c()
V <- c()
W <- c()
X <- c()
Y <- c()
stop <- c()

for(i in accs) {
  #get genbank data using accession number as i
  id = GBAccession(i)

  #read in genbank sequence data, using partial as true to include all cds, not just complete ones
  gb = readGenBank(id, partial = TRUE)
  #subset genbank data to only the gene identified as tetM
  gb_cds_tetM <- subset(gb@cds, gene=="tetM")
  
  #identify start and end location of coding strand
  start <- gb_cds_tetM@ranges@start
  end <- gb_cds_tetM@ranges@width+gb_cds_tetM@ranges@start-1
  #get DNAString for sequence
  seq <- gb@sequence
  #subset the whole sequence with the coding strand that codes for tetM
  cds_seq <- subseq(seq, start = start, end = end)
  
  #create list of accession numbers in conjuction with gc content
  Accession <- c(Accession, i)
  #calculate GC percentage
  GC <- c(GC, ((alphabetFrequency(cds_seq)[, 2]/width(cds_seq))+(alphabetFrequency(cds_seq)[, 3]/width(cds_seq))))

  aa_seq <- gb_cds_tetM$translation
  aa_freq <- c(aa_freq, alphabetFrequency(gb_cds_tetM$translation))
}

GCcont <- data.frame(Accession, GC)

