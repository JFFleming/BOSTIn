#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
LB <- read.table(filename, sep="\t", header=TRUE)

#Make the histogram
summ_pdf <- paste0(tools::file_path_sans_ext(filename),".histogram.pdf")
pdf(file = summ_pdf, width=4, height=4)
hist(LB[,2], main="Histogram of Genewise LB Score values", xlab="LB Score Values")
dev.off()

#Calculate stats
summ_mean = mean(LB[,2])
summ_median = median(LB[,2])
summ_sd = sd(LB[,2])
summ_mad = mad(LB[,2])
summ_min = min(LB[,2])
summ_max = max(LB[,2])
summ_smallmad = summ_median+(2*summ_mad)
summ_bigmad = summ_median+(3*summ_mad)


#Make the Red and Yellow Lists
redfile <- paste0(tools::file_path_sans_ext(filename),".BranchLength.GeneWise.redflags.txt")
summ_RedFlags<-file(redfile)
redList <- c()
yellowfile <- paste0(tools::file_path_sans_ext(filename),".BranchLength.GeneWise.yellowflags.txt")
summ_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(LB[,2])){
	if(LB[i,2] >= summ_bigmad){
		redGene <- paste(LB[i,1], "\t", LB[i,2], "\n")
		redList<- append(redList, redGene)
		yellowGene <- paste(LB[i,1], "\t", LB[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}
	else if(LB[i,2] >= summ_smallmad){
		yellowGene <- paste(LB[i,1], "\t", LB[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}

}
cat(redList, file=summ_RedFlags)
cat(yellowList, file=summ_YellowFlags)
#Write the Narrative Summary
summ_file <- paste0(tools::file_path_sans_ext(filename),".NarrativeSummary.txt")
summ_list <- c()

summ_list <- append(summ_list, "Welcome to the GeneWise batch selection narrative summary for the LB Score, the branch length heterogeneity metric used by BOSTIn.\n")
summ_list <- append(summ_list, paste("You can find a summary histogram of the standard deviation of the LB Score for each gene here:\n",summ_pdf,"\n"))
summ_list <- append(summ_list, "\n\nGenewise LB Score summary statistics\n")

summ_list <- append(summ_list, paste("mean\t", summ_mean, "\n"))
summ_list <- append(summ_list, paste("median\t", summ_median, "\n"))
summ_list <- append(summ_list, paste("standard deviation\t", summ_sd, "\n"))
summ_list <- append(summ_list, paste("median absolute deviation\t", summ_mad, "\n"))
summ_list <- append(summ_list, paste("minimum\t", summ_min, "\n"))
summ_list <- append(summ_list, paste("maximum\t", summ_max, "\n"))
summ_list <- append(summ_list, "\n\nThe LB-Score measures the percentage deviation of each taxon from the average patristic distance, and so is independent of the root of the tree, making it quite useful to identify long branches.\nBostIn rapidly generates a Neighbour-Joining tree to calculate the LB-Score. This produces an LB-Score that is normally significantly similar, even under large amounts of Long Branch Attraction, but it won't be as accurate as an LB-Score generated under the best possible model.\nFor the purposes of defining the sextile of genes most likely to cause a long branch attraction artifact, however, it ought to suffice.\n")
summ_list <- append(summ_list, paste("The LB-Scores in your dataset range from ",summ_min," to ",summ_max,", with a mean of",summ_mean,"and a standard deviation of",summ_sd,".\n"))
summ_list <- append(summ_list, "Using the median and the median absolute deviation, which is similar to a standard deviation, but for medians, we can use the LB-Score to more robustly identify suspect genes with large amounts of branch length heterogeneity by assessing which genes are outside of two, and then 3 median absolute deviations of the median of all the genes in the dataset.\nThis is because it is the extremes of branch length heterogeneity that can cause the greatest problems.")
summ_list <- append(summ_list, paste("The median absolute deviation is",summ_median,"and so the median plus two median absolute deviations is",summ_smallmad,", while plus three median absolute deviations is ",summ_bigmad,".\nWe've identified genes beyond these bounds as Yellow Flags and Red Flags respectively, as with the other measurements in BOSTIn.\nYour Red Flag genes are:\n"))
summ_list <- append(summ_list, paste(redList, sep="\n"))
summ_list <- append(summ_list, "\nYour Yellow Flag genes are:\n")
summ_list <- append(summ_list, paste(yellowList, sep="\n"))
cat(summ_list, file=summ_file)
