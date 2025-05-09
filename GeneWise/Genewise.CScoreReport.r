#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
SAT <- read.table(filename, sep="\t", header=TRUE)

#Make the histogram
summ_pdf <- paste0(tools::file_path_sans_ext(filename),".histogram.pdf")
pdf(file = summ_pdf, width=4, height=4)
hist(SAT[,2], main="Histogram of Genewise C Score values", xlab="C Score Values")
dev.off()

#Calculate stats
summ_mean = mean(SAT[,2])
summ_median = median(SAT[,2])
summ_sd = sd(SAT[,2])
summ_mad = mad(SAT[,2])
summ_min = min(SAT[,2])
summ_max = max(SAT[,2])
summ_smallmad = summ_median+(2*summ_mad)
summ_bigmad = summ_median+(3*summ_mad)


#Make the Red and Yellow Lists
redfile <- paste0(tools::file_path_sans_ext(filename),".SiteSaturation.GeneWise.redflags.txt")
summ_RedFlags<-file(redfile)
redList <- c()
yellowfile <- paste0(tools::file_path_sans_ext(filename),".SiteSaturation.GeneWise.yellowflags.txt")
summ_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(SAT[,2])){
	if(SAT[i,2] <= 10){
		redGene <- paste(SAT[i,1], "\t", SAT[i,2], "\n")
		redList<- append(redList, redGene)
		yellowGene <- paste(SAT[i,1], "\t", SAT[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}
	else if(SAT[i,2] <= 20){
		yellowGene <- paste(SAT[i,1], "\t", SAT[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}

}
cat(redList, file=summ_RedFlags)
cat(yellowList, file=summ_YellowFlags)
#Write the Narrative Summary
summ_file <- paste0(tools::file_path_sans_ext(filename),".SiteSaturation.NarrativeSummary.txt")
summ_list <- c()

summ_list <- append(summ_list, "Welcome to the GeneWise batch selection narrative summary for the C Score, the nucleotide saturation metric used by BOSTIn.\n")
summ_list <- append(summ_list, paste("You can find a summary histogram of the standard deviation of the C Score for each gene here:\n",summ_pdf,"\n"))
summ_list <- append(summ_list, "\n\nGenewise C Score summary statistics\n")

summ_list <- append(summ_list, paste("mean\t", summ_mean, "\n"))
summ_list <- append(summ_list, paste("median\t", summ_median, "\n"))
summ_list <- append(summ_list, paste("standard deviation\t", summ_sd, "\n"))
summ_list <- append(summ_list, paste("median absolute deviation\t", summ_mad, "\n"))
summ_list <- append(summ_list, paste("minimum\t", summ_min, "\n"))
summ_list <- append(summ_list, paste("maximum\t", summ_max, "\n"))
summ_list <- append(summ_list, "\n\nAs you are using nucleotide data, we selected the C Factor - the Convergence Factor. This is based on the ratio of the standard deviation of the Transition/Transversion ratio of the dataset and the standard deviation of the uncorrected p-distance.\n")
summ_list <- append(summ_list, paste("The C Scores in your dataset range from ",summ_min," to ",summ_max,", with a mean of",summ_mean,"and a standard deviation of",summ_sd,".\n"))
summ_list <- append(summ_list, "We can identify genes subjective to significant amounts of saturation using the C Score relatively easily. Genes with C Scores below 10 are marked as Red Flag genes, while those with C Scores below 20 are marked as Yellow Flag genes. Your Red Flag genes are:\n")
summ_list <- append(summ_list, paste(redList, sep="\n"))
summ_list <- append(summ_list, "\nYour Yellow Flag genes are:\n")
summ_list <- append(summ_list, paste(yellowList, sep="\n"))
cat(summ_list, file=summ_file)
