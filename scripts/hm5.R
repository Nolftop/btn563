library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(scales)


dna <- fasta2DNAbin(file="/home/ns/Documents/classes/botany/btn563/results/align/out.fasta")

D <- dist.dna(dna, model="TN93")

M <-as.matrix(D)

colnames(M) <- strtrim(colnames(M), 10)
rownames(M) <- strtrim(rownames(M), 10)
maxs <- apply(M, 2, max)
mins <- apply(M, 2, min)
M <- rescale(M)
M[is.nan(M)] <- 0

temp <- as.data.frame(M)
temp

plot <- table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)

