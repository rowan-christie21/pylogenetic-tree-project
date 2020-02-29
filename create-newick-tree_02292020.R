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
library(rgl)

#define path to sequence
path <- "C:/Users/smc9227/Documents/bioinformatics_code_stuff/antibiotic_resistant_strains3.fasta"

#_____________________
#1st Step: Perform MSA

#Function for performing msa
msa_data <- function(path_to_seq=path, format = "fasta") {
  #call to compiled fasta file with all of the bacterial sequences
  seq <- readDNAStringSet(path_to_seq, format = format)
  #call MSA function which performs multiple sequence alignment
  alignment <- msa(seq)
  return(alignment)
}

#works!
msa_data()

#_____________________

#2nd Step: Convert the MSA to PhyDat class so it can be used in a tree

convert_msa_to_phyDat <- function(alignment=msa_data()) {
  msa_phyDat <- msaConvert(alignment, type="phangorn::phyDat")
  return(msa_phyDat)
}

#works!
convert_msa_to_phyDat()

#_________________________________
#3rd Step: Create distance matrix

#Default model is Jukes Cantor model (1969), which is the 
#simplest model of evolution assumes that each base has an equal chance of changing
distance_matrix <- function(phyDat=convert_msa_to_phyDat(), model="JC69") {
  dm <- dist.ml(phyDat, model=model)
  return(dm)
}

#works!
distance_matrix()

#__________________________________________
#4th Step: Estimate trees using algorithms

#estimation algorithims from phangorn:
#https://cran.r-project.org/web/packages/phangorn/phangorn.pdf
  #upgma                   UPGMA and WPGMA  
  #NJ                      Neighbor-Joining
  #neighborNet             Computes a neighborNet from a distance matrix
  
estimation_algo <- function(dm=distance_matrix()) {
  upgma_algo <- upgma(dm)
  nj_algo <- nj(dm)
  return(list(upgma_algo, nj_algo))
}

estimation_algo.plot <- function(est_algo=estimation_algo()) {
  plot(est_algo[[1]], main="UPGMA")
  plot(est_algo[[2]], main="NJ")
}

#works!
estimation_algo()

#works!
estimation_algo.plot()

#______________________________________
#5th Step: Determine maximum parsimony

max_parsimony <- function(est_algo=estimation_algo(), phyDat=convert_msa_to_phyDat()) {
  phylo_upgma <- est_algo[[1]]
  phylo_nj <- est_algo[[2]]
  
  upgma_pars <- parsimony(phylo_upgma, phyDat)
  nj_pars <- parsimony(phylo_nj, phyDat)
  
  pars_data <- data.frame(c("UPGMA", upgma_pars), c("Neighbor-Joining", nj_pars))

  if(upgma_pars > nj_pars) {
    print(pars_data)
    print("UPGMA has maximum parsimony.")
  }
  if(nj_pars > upgma_pars) {
    print(pars_data)
    print("Neighbor-Joining has maximum parsimony.")
  }
  
  phylo.input <- readline(prompt = "Which model would you like to use: UPGMA or NJ? ")
  
  if(phylo.input == "UPGMA") {
    return(phylo_upgma)
  }
  if(phylo.input == "NJ") {
    return(phylo_nj)
  }
  
  #return(pars_data)
}

pars <- max_parsimony(upgma_test, nj_test, phyDat_test)


#______________________________________
#6th Step: Maximum likelihood approach

#default model is Jukes Cantor
max_likelihood <- function(phylo, phyDat, model="JC", rearrangement="stochastic") {
  #computing the likelihood of a given tree
  fit <- pml(phylo, phyDat)
  print(fit)
  
  #optimize tree topology and branch length for a given model
  optimize <- optim.pml(fit, model= model, rearrangement = rearrangement)
  return(optimize)
}

ml_test <- max_likelihood(pars, phyDat_test)


#_____________________

#master function for creating newick tree
create_tree <- function(path, model) {
  msa <- msa_data(path)
  phyDat <- convert_msa_to_phyDat(msa)
  dm <- distance_matrix(phyDat)
  
  #right now we only have phyDat function so far so return that
  return(dm)
}

path <- "C:/Users/smc9227/Documents/bioinformatics_code_stuff/antibiotic_resistant_strains3.fasta"
tree_test <- create_tree(path)

















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


