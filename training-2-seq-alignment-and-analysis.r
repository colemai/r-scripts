#!/usr/bin/Rscript
# Purpose: Seq alignment and Analysis
# Author: Ian Coleman
#rm(list = ls())

# Fetch bioconductor packages
source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

# Libraries
library(Biostrings)
library(seqinr)

## ================== Alignment and Save to FASTA ========================= ##

# Read in the fasta
prokaryotes <- read.fasta(file="fastas/prok.fasta",seqtype="DNA")

# Pull out the first two seqs as characters
seq1 <- as.character(prokaryotes[[1]])
seq1 <- paste(seq1, collapse="")

seq2 <- as.character(prokaryotes[[2]])
seq2 <- paste(seq2, collapse="")

# Pairwise alignment with the default function
pairAlign <- pairwiseAlignment(pattern=seq2, subject=seq1)
summary(pairAlign)

# Write this alignment to a new FASTA file
# First convert to string
pairAlignString <- BStringSet(c(toString(subject(pairAlign)), toString(pattern(pairAlign))))
writeXStringSet(pairAlignString, "aligned.txt", format="FASTA")


## ================== DotPlot Alignment ========================= ##

# Get two protein sequences (less noisy on dot plot than NTs)
coxGenes <- read.fasta(file="fastas/cox1multi.fasta", seqtype="AA")
cox1 <- as.character(coxGenes[[1]])
cox2 <- as.character(coxGenes[[2]])

dotPlot(cox1, cox2, main="Human vs Mouse COX1 Dotplot")

# a windowed dotplot is easier to interpret:
dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3,
        wstep = 3, nmatch = 3")

# First 100 residues only
dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3,
        wstep = 3, nmatch = 3")

# Three main types of alignment (i) Global (ii) Local (iii) Overlap - for overlapping seqs
