#3/12/2020
#trying to learn how to retrieve sequences and corresponding
#sequence features (like cds, protein translation) from genbank

library(rvest)    
scotusURL <- "https://en.wikipedia.org/wiki/Novosibirsk"
  
## ********************
## Option 1: Grab the tables from the page and use the html_table function to extract the tables you're interested in.

temp <- scotusURL %>% 
  html %>%
  html_nodes("table")

html_table(temp[1]) ## Just the "legend" table
html_table(temp[2]) ## THE MAIN TABLE


library(AnnotationBustR)
library(reutils)

vignette("AnnotationBustR-vignette")

test.search <- esearch(term = "tetM", db="nuccore") #nuccore = nucleotide

accessions<-efetch(test.search, rettype = "acc",retmode = "text")#fetch accessions

accessions <- strsplit(content(accessions), "\n")[[1]]#split out accessions from other meta-data

#use =textjoin(""", """, 1, E6:E34) to copy cells
accessions <- c("MG874353", "MG874352", "MG874351", "MG874350", "MG874349", "MG874348", "MG874347", "MG874346", "MG874345", "MG874344", "MG874342", "MG874341", "MG874340", "MG874339", "MG874338", "MG874337", "MG874336", "MG874335", "MG874334", "MG874333", "MG874332", "MG874331", "MG874330", "MG874329", "MG874328", "MG874327", "MG874326", "MG874325", "MG874343")
  
acc.search <- esearch(term = accessions, db="nuccore")

Locus <- c(accessions)
term_df <- data.frame(Locus)
term_df$Type <- "CDS"
term_df$Name <- "tetM"


install.packages("ape")
install.packages("seqinr")

library(ape) #this is a general R-package for phylogenetics and comparative methods
library("seqinr") #this is an specialized package for nucleotide sequence management

n.gon <- c("MG874353", "MG874352", "MG874351", "MG874350", "MG874349", "MG874348", "MG874347", "MG874346", "MG874345", "MG874344", "MG874342", "MG874341", "MG874340", "MG874339", "MG874338", "MG874337", "MG874336", "MG874335", "MG874334", "MG874333", "MG874332", "MG874331", "MG874330", "MG874329", "MG874328", "MG874327", "MG874326", "MG874325", "MG874343")

n.gon.seq <- read.GenBank(n.gon)

str(n.gon.seq)
attributes(n.gon.seq)
attr(n.gon.seq, "species")


seqs <- AnnotationBust(Accessions = accessions, Terms = term_df)


library (rentrez)

n_search <- entrez_search(db="nuccore", term="MG874353") 
n_search$ids

protein_links <- entrez_link(dbfrom='nuccore', id=n_search$ids, db='all')
protein_seq <- entrez_fetch(db="protein", rettype="fasta", id=protein_links$links$nuccore_protein)


library(UniprotR)

n_id <- UniProt.ws(taxId=9606)

n.gon <- c("MG874353", "MG874352", "MG874351", "MG874350", "MG874349", "MG874348", "MG874347", "MG874346", "MG874345", "MG874344", "MG874342", "MG874341", "MG874340", "MG874339", "MG874338", "MG874337", "MG874336", "MG874335", "MG874334", "MG874333", "MG874332", "MG874331", "MG874330", "MG874329", "MG874328", "MG874327", "MG874326", "MG874325", "MG874343")

n.gon.taxa <- GetNamesTaxa(n.gon,"list.csv")


n_search <- entrez_search(db="nuccore", term="MG874353") 
n_fetch <- entrez_fetch(db = "nuccore", id = n_search$ids, rettype = "gb")


n_links <- entrez_link(dbfrom='nuccore', id=n_search$ids, db='all')

protein_sum <- entrez_summary(db="protein", rettype="fasta", id=n_links$links$nuccore_protein)

protein_seq <- entrez_fetch(db="protein", rettype="gb", id=n_links$links$nuccore_protein)
protein_seq_for <- strsplit(protein_seq, ">")[[1]][1]

library(biomaRt)
library(biomartr)


library(read.gb)

#read.gb("https://www.ncbi.nlm.nih.gov/nuccore/MG874353.1")
n_gb <- read.gb("n_gon.gb", Type = "full")
n_gb$MG874353$FEATURES

#----------------------------------------------------------
#will use genbankr to retreive sequence info

library("genbankr")

n_genbank <- readGenBank("n_gon.gb")


id1 = GBAccession("JF915701")
gb1 = readGenBank(id1)
gb1@cds$gene
gb1_cds_tetM <- subset(gb1@cds, gene=="tetM")
start1 <- gb1_cds_tetM@ranges@start
end1 <- gb1_cds_tetM@ranges@width+gb1_cds_tetM@ranges@start-1
seq1 <- gb1@sequence
cds_seq1 <- subseq(seq1, start = start1, end = end1)

gb1_cds_tetM$translation
cds_aa_seq1 <- Biostrings::translate(cds_seq1)

id2 = GBAccession("X92947")
gb2 = readGenBank(id2)
gb2@cds$gene
gb2_cds_tetM <- subset(gb2@cds, gene=="tetM")
start2 <- gb2_cds_tetM@ranges@start
end2 <- gb2_cds_tetM@ranges@width+gb2_cds_tetM@ranges@start-1
seq2 <- gb2@sequence
cds_seq2 <- subseq(seq2, start = start2, end = end2)

gb2_cds_tetM$translation
cds_aa_seq2 <- Biostrings::translate(cds_seq2)

#get genbank info from accession number, partial is true to allow
#for partial cds
id3 = GBAccession("MG874333")
gb3 = readGenBank(id3, partial = TRUE)

#subset genes/products by tetM term
gb3@cds$gene
gb3_cds_tetM <- subset(gb3@cds, gene=="tetM")

#get start and end of coding sequence (cds)
start3 <- gb3_cds_tetM@ranges@start
end3 <- gb3_cds_tetM@ranges@width+gb3_cds_tetM@ranges@start-1

#subset cDNA sequence by coding sequence range (start and end)
seq3 <- gb3@sequence
cds_seq3 <- subseq(seq3, start = start3, end = end3)

#double check everything is subsetted correctly and results
#translate same product
gb3_cds_tetM$translation
cds_aa_seq3 <- Biostrings::translate(cds_seq3)

seqtest <- c(cds_seq1, cds_seq2, cds_seq3)

msa_res <- msa(seqtest)
print(msa_res)
msaPrettyPrint(msa_res)
texi2pdf("msa_res.tex", clean=TRUE)

data(BLOSUM62)
scores <- msaConservationScore(msa_res, BLOSUM62)

msa_mask <- msa_res
colM <- IRanges(start = 1, end = 100)
colmask(msa_mask) <- colM

unmasked(msa_mask)


conmat <- consensusMatrix(msa_res)
conmat[, 101:110]

