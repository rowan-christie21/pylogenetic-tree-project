# Feb 29 2020
# File including functions for creating a phylogentic tree in
# Newick format
# Defaults are:

#_____________________
# Load Libraries

library(ape)
library(phangorn)
library(seqinr)
library(BiocManager)
library(Biostrings)
library(msa)
library(ggplot2)
library(ggtree)

#_____________________
#1st Step: Perform MSA

#Function for performing msa
msa_data <- function(path_to_seq, format = "fasta") {
  #call to compiled fasta file with all of the bacterial sequences
  seq <- readDNAStringSet(path_to_seq, format = format)
  #call MSA function which performs multiple sequence alignment
  alignment <- msa(seq)
  return(alignment)
}

#works!
msa_test <- msa_data("C:/Users/smc9227/Documents/bioinformatics_code_stuff/antibiotic_resistant_strains3.fasta", "fasta")

#_____________________

#2nd Step: Convert the MSA to PhyDat class so it can be used in a tree

convert_msa_to_phyDat <- function(alignment) {
  msa_phyDat <- msaConvert(alignment, type="phangorn::phyDat")
  return(msa_phyDat)
}

phyDat_test <- convert_msa_to_phyDat(msa_test)



#master function for creating newick tree
create_phyDat <- function(path, model) {
  the_msa <- msa_data(path)
  the_phyDat <- convert_msa_to_phyDat(the_msa)
  
  #right now we only have phyDat function so far so return that
  return(the_phyDat)
}

path <- "C:/Users/smc9227/Documents/bioinformatics_code_stuff/antibiotic_resistant_strains3.fasta"

create_phyDat(path)




#perform conversion function to phyDat class so phangorn can operate on it
bact_phyDat <- msaConvert(bact_alignment, type="phangorn::phyDat")
















#---------------------------------------------------------

# PYLOGENETIC TREE EXAMPLE DEC 19 2020

#---------------------------------------------------------
#Process to create a phylogenetic tree for
#bacterial genes with TetM gene
#Dec 14 2020

library(ape)
library(phangorn)
library(seqinr)
library(BiocManager)
library(Biostrings)
library(msa)
library(ggplot2)
library(ggtree)

#1st Step: Perform MSA
#_____________________

#call to compiled fasta file with all of the bacterial sequences
bact_seq <- readDNAStringSet("C:/Users/smc9227/Documents/bioinformatics_code_stuff/antibiotic_resistant_strains3.fasta", format = "fasta")

#call MSA function which performs multiple sequence alignment
bact_alignment <- msa(bact_seq)

#check results
print(bact_alignment, show="complete")


#2nd Step: Convert the MSA to PhyDat class so it can be used in a tree
#_____________________________________________________________________

#perform conversion function to phyDat class so phangorn can operate on it
bact_phyDat <- msaConvert(bact_alignment, type="phangorn::phyDat")

#3rd Step: Create distance matrix
#_________________________________

#use Jukes Cantor model (1969), which is the simplest model of evolution
#assumes that each base has an equal chance of changing
dm <- dist.ml(bact_phyDat, model="JC69")

#4th Step: Estimate trees using algorithms
#__________________________________________

#using neighbor-joining and UPGMA algorithms
#to estimate trees from distance matrices
bact_UPGMA <- upgma(dm)
bact_NJ  <- NJ(dm)

#visualize trees
plot(bact_UPGMA, main="UPGMA")
plot(bact_NJ, main = "Neighbor Joining")


#5th Step: Determine maximum parsimony
#______________________________________

#find out which tree is most parsimonous to find out 
#which model is the better fit for the data
parsimony(bact_UPGMA, bact_phyDat)
#parsimony score: 1960
parsimony(bact_NJ, bact_phyDat)
#parsimony score: 1943
#bact_UPGMA has maximum parsimony, so we will use that

#6th Step: Maximum likelihood approach
#______________________________________

#computing the likelihood of a given tree
fit <- pml(bact_UPGMA, bact_phyDat)
print(fit)
#optimize tree topology and branch length for Jukes Cantor model
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")

#7th step: Visualize tree
#________________________

#Using ggplot's package ggtree to visualize it
ggtree(fitJC$tree) +
  geom_tiplab(hjust = -.04) +
  xlim(0, 1) +
  geom_tippoint(size=2, color="blue")

#---------------------------------------------
#Now try to write tree to file
write.tree(fitJC$tree, "fitJC.txt")


