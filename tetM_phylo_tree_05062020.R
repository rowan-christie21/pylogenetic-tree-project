#5/6/2020
#Creating phylogenetic tree of bacteria

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
library(ape)
library(phangorn)
library(seqinr)
library(BiocManager)
library(Biostrings)
library(msa)
library(ggtree)


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
  
  #add sequences to vector in character form to organize it into a list of dnastringsets more easily
  seqs <- c(seqs, as.character(cds_seq))
}

#convert seqs vector into list of dnastringsets
seqlist <- DNAStringSet(paste0(as.character(seqs)))

#perform multiple sequence alignment with all sequences
tetm_tree <- msa(seqs, type = "dna")
#print multiple sequence alignment
print(tetm_tree)
#print to pdf
msaPrettyPrint(tetm_tree)
#if pdflatex does not work, use TeXworks app to print tex file

#perform conversion function to phyDat class so phangorn can operate on it
tetm_phydat <- msaConvert(tetm_tree, type="phangorn::phyDat")

#use Jukes Cantor model (1969), which is the simplest model of evolution
#assumes that each base has an equal chance of changing
dm <- dist.ml(tetm_phydat, model="JC69")

#using neighbor-joining and UPGMA algorithms
#to estimate trees from distance matrices
tetm_UPGMA <- upgma(dm)
tetm_NJ  <- NJ(dm)

#visualize trees
plot(tetm_UPGMA, main="UPGMA")
plot(tetm_NJ, main = "Neighbor Joining")

#find out which tree is most parsimonous to find out 
#which model is the better fit for the data
parsimony(tetm_UPGMA, tetm_phydat)
#parsimony score: 329
parsimony(tetm_NJ, tetm_phydat)
#parsimony score: 318
#tetm_UPGMA has maximum parsimony, so we will use that

#computing the likelihood of a given tree
fit <- pml(tetm_UPGMA, tetm_phydat)
print(fit)
#optimize tree topology and branch length for Jukes Cantor model
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")

#add missing strain name
fitJC$tree$tip.label[1] <- "Kagoshima6"

#add dataframe for labeling tree
dd <- data.frame(strain=tetm_strains$Strain, Species=tetm_strains$Species, MIC=tetm_strains$`tetracycline µg/ml`)

#plot tree
p <- ggtree(fitJC$tree)

#annotate tree
p <- p %<+% dd + geom_tiplab(aes(color=Species)) + 
  geom_tippoint(aes(size=MIC, color=Species), alpha=0.25) +
  xlim(0,0.15) +
  theme(legend.position="right") + 
  labs(size="Tetracycline µg/ml") +
  geom_strip('Kaluga 10219', 'Kazan 09/15/001', barsize=2, color='purple', label="Clade 1", offset.text=.005, offset=0.03) +
  geom_strip('CHS-1E', 'JH2-2', barsize=2, color='red', label="Clade 2", offset.text=.005, offset=0.03) +
  geom_strip('Bryansk 36/16/08', 'Novosibirsk 10226', barsize=2, color='blue', label="Clade 3", offset.text=.005, offset=0.03) +
  geom_strip('RV-16', 'DMV42A', barsize=2, color='green', label="Clade 4", offset.text=.005, offset=0.03)

p

png(filename = paste("D:/Bioinformatics-Projects/phylogentic-tree-project-02292020/tetm_tree ",Sys.Date(),".png", sep = ''), width = 900, height = 600)
p
dev.off()



