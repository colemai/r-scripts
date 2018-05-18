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

## ================== FASTA Processing - NT characteristics ========================= ##
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
n <- length(GCwindow) # how many windows are there
Chunks <- numeric(n) # make an empty variable of that size

#For loop to compute GC content per chunk
for (i in 1:n) {
  chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
  chunkGC <- GC(chunk)
  Chunks[i] <- chunkGC
}

## ================== Plotting the FASTA characteristics ========================= ##

plot(GCwindow, Chunks, type="b", xlab="NT start position", ylab="GC content")


## ================== Create a function to display GC content ========================= ##
# function to take seq and window size and plot GC content
slidingwindowGCplot <- function(windowsize, inputseq)
{
  GCwindow <- seq(1, length(inputseq)-windowsize, by=windowsize)
  n <- length(GCwindow) # find window size
  Chunks <- numeric(n) # make blank vector of that size
  
  #iterate through the windows of the seq and calc GC content
  for (i in 1:n) {
    chunk <- inputseq[GCwindow[i]:(GCwindow[i]+windowsize-1)]
    Chunks[i] <- GC(chunk)
  }
  
  #plot the GC content of the windows
  plot(GCwindow, Chunks, type="b", xlab="NT start pos", ylab = "GC Content")
  
}

