#!/usr/local/bin/Rscript

library(phangorn)
args <- commandArgs(trailingOnly = TRUE)
a_file <- args[1]
a_tree <- args[2]
at_out <- args[3]
a_out <-args[4]
d_type <- args[5]
Fasta=read.phyDat(a_file, format = "fasta", type = d_type)
DistMatrix <- dist.ml(Fasta)
treeNJ <- NJ(DistMatrix)
write.tree(treeNJ, file=a_tree)

dist.mat = cophenetic.phylo(treeNJ)
dist.mat[dist.mat == 0] <- NA
matrixmean <- mean(dist.mat, na.rm=TRUE)
result_matrix <- matrix(NA, nrow = nrow(dist.mat), ncol = 2)

for (col in 1:ncol(dist.mat)) {
   t_avgdist <- mean(dist.mat[col,], na.rm=TRUE)
   t_deviation <- ((t_avgdist/matrixmean)-1)*100
   result_matrix[col, ] <- c(t_avgdist, t_deviation)
}

rownames(result_matrix) <- rownames(dist.mat)
colnames(result_matrix) <- c("averagePDist", "LB-Score")
UQ <- result_matrix[,2][result_matrix[,2] >= quantile(result_matrix[,2])[4]]
UQ

write.table(result_matrix, sep="\t", row.names=TRUE, file=at_out)
write(paste0("standardDeviation\t", sd(result_matrix[,2])), file=a_out)
write(paste0("Mean\t", mean(result_matrix[,2])), file=a_out, append=TRUE)
write(paste0("Minimum\t", quantile(result_matrix[,2])[1]), file=a_out, append=TRUE)
write(paste0("LowerQuartile\t", quantile(result_matrix[,2])[2]), file=a_out, append=TRUE)
write(paste0("Median\t", quantile(result_matrix[,2])[3]), file=a_out, append=TRUE)
write(paste0("UpperQuartile\t", quantile(result_matrix[,2])[4]), file=a_out, append=TRUE)
write(paste0("Maximum\t", quantile(result_matrix[,2])[5]), file=a_out, append=TRUE)
write(paste0("stdDevUpperQuartile\t", sd(UQ)), file=a_out, append=TRUE)

t_pdf <- paste0(tools::file_path_sans_ext(a_out),".histogram.pdf")
pdf(file = t_pdf, width=4, height=4)
#hist(result_matrix[,2], main="Histogram of LB-Scores", xlab="LB-Scores")
hist(UQ, main="Histogram of Upper Quartile LB-Scores", xlab="LB-Scores")
dev.off()