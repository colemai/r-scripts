#!/usr/bin/Rscript
# An R Script to work with Fasta files, to learn how to do so
# Ian Coleman

## ================== Import Local Fasta ========================= ##
#install.packages("seqinr")
library(seqinr)

cox1 <- read.fasta(file="fastas/cox1.fasta", seqtype="AA")
length(cox1) #4
length(cox1[1])

seq1 <- cox1[1]

## ================== Import GenBank Fasta ========================= ##
#install.packages("ape")
library(ape)

# Retrieve a GenBank seq as a binary object
AB003468 <- read.GenBank("AB003468", as.character = "TRUE") 

# Save Genbank Sequence as FASTA format
write.dna(AB003468, file="AB003648.fasta", format="fasta", append=
            FALSE, nbcol=6, colsep="", colw=10)

## ================== Import GenBank Fasta via rentrez package ========================= ##
# Rentrez allows us to accress data in formats other than binary
# and by means other than accession number
install.packages("rentrez") 
library(rentrez)

# Search for a seq, this returns many hits tho so sticking with ape
entrez_search(db="nucleotide", term="human superoxide dismutase")

## ================== FASTA Processing ========================= ##
#Take the first sequence (there's only one but this removes header)
CloningVector <- AB003468[[1]]

#Count each of a,c,g,t incidence
count <- count(CloningVector, 1)

#Count each 2NT combination
count2 <- count(CloningVector, 2)

#Get the GC Content (G+C:A+T) - as coding regions generally higher in GC {seqinr lib}
GC <- GC(CloningVector)

#That was only the global GC though, let's look at smaller patches. Window of 200
GCwindow <- seq(1,length(CloningVector)-200, by=200)
n <- length(GCwindow)
Chunks <- numeric(n)

#For loop to compute GC content per chunk
for (i in 1:n) {
  chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
  chunkGC <- GC(chunk)
  Chunks[i] <- chunkGC
}
