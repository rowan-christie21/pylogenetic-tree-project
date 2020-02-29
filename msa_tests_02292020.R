# Nov 8
# testing out bioconductor package msa

library(msa)

system.file("tex", "texshade.sty", package="msa")

mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences

myFirstAlignment <- msa(mySequences)

## use default substitution matrix
myFirstAlignment

print(myFirstAlignment, show="complete")

msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

# restrict to a teasing example
msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis",
               showNames="none", showLogo="none", askForOverwrite=FALSE)

#muscle alignment
muscle_alignment <- msa(mySequences, "Muscle")

## use default substitution matrix
muscle_alignment

msaPrettyPrint(muscle_alignment, output="pdf",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

myMaskedAlignment <- muscle_alignment
colM <- IRanges(start=1, end=100)
colmask(myMaskedAlignment) <- colM
myMaskedAlignment

unmasked(myMaskedAlignment)

conMat <- consensusMatrix(myFirstAlignment)
dim(conMat)


####
#test conversion

bact_seq <- readDNAStringSet("C:/Users/smc9227/Documents/bioinformatics_code_stuff/antibiotic_resistant_strains2.fasta", format = "fasta")

bact_alignment <- msa(bact_seq)

print(bact_alignment, show="complete")

msaPrettyPrint(bact_alignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

#converting to phyDat class
bact_phyDat <- msaConvert(bact_alignment, type="phangorn::phyDat")

mt <- modelTest(bact_phyDat)
print(mt)
dm <- dist.ml(bact_phyDat, model="JC69")

mammals_UPGMA <- upgma(dm)

#plot(mammals_UPGMA, main="UPGMA")

#making tree prettier
library("ggplot2")
library("ggtree")

#tree <- read.tree(nwk)
mammals_UPGMA$tip.label
groupInfo <- split(mammals_UPGMA$tip.label, gsub("Neisseria", "N.", gsub(", partial cds", "", mammals_UPGMA$tip.label)))
mammals_UPGMA <- groupOTU(mammals_UPGMA, groupInfo)

#mammals_UPGMA[3]$tip.label

#r <- c(1:23)
#ggtree(mammals_UPGMA)$data$r

ggtree(mammals_UPGMA) + #+ geom_tiplab() + geom_treescale() + xlim(0, 0.06)
  geom_tiplab(hjust = -.01) +
  xlim(0, 1)
  #geom_tiplab(hjust = -.01) + xlim(0, .05)
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.5)
  #geom_point(aes(fill = rate), shape = 21, size = 4) +
  #scale_color_manual(values = c("black", "red"), guide = FALSE) +
  #scale_fill_continuous(low = 'blue', high = 'red') +
  #theme_tree2()# + theme(legend.position = 'right')

#--

#converting to seqinr class
hemoAln2 <- msaConvert(bact_alignment, type="seqinr::alignment")

library(seqinr)
d <- dist.alignment(hemoAln2, "identity")
#as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

library(ape)
hemoTree <- nj(d)
plot(hemoTree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")

####

#test bootstrap
write.tree(bs, file="bootstrap_example.tre")


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


#---


