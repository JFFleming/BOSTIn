#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
CH <- read.table(filename, sep="\t", header=TRUE)

#Make the histogram
summ_pdf <- paste0(tools::file_path_sans_ext(filename),".histogram.pdf")
pdf(file = summ_pdf, width=4, height=4)
hist(CH[,2], main="Histogram of Genewise nRCFV values", xlab="nRCFV Values")
dev.off()

#Calculate stats
summ_mean = mean(CH[,2])
summ_median = median(CH[,2])
summ_sd = sd(CH[,2])
summ_mad = mad(CH[,2])
summ_min = min(CH[,2])
summ_max = max(CH[,2])
summ_smallmad = summ_median+(2*summ_mad)
summ_bigmad = summ_median+(3*summ_mad)


#Make the Red and Yellow Lists
redfile <- paste0(tools::file_path_sans_ext(filename),".CompHet.GeneWise.redflags.txt")
summ_RedFlags<-file(redfile)
redList <- c()
yellowfile <- paste0(tools::file_path_sans_ext(filename),".CompHet.GeneWise.yellowflags.txt")
summ_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(CH[,2])){
	if(CH[i,2] >= summ_bigmad){
		redGene <- paste(CH[i,1], "\t", CH[i,2], "\n")
		redList<- append(redList, redGene)
		yellowGene <- paste(CH[i,1], "\t", CH[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}
	else if(CH[i,2] >= summ_smallmad){
		yellowGene <- paste(CH[i,1], "\t", CH[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}

}
cat(redList, file=summ_RedFlags)
cat(yellowList, file=summ_YellowFlags)
#Write the Narrative Summary
summ_file <- paste0(tools::file_path_sans_ext(filename),".CompHet.NarrativeSummary.txt")
summ_list <- c()

summ_list <- append(summ_list, "Welcome to the GeneWise batch selection narrative summary for the nRCFV, the compositional heterogeneity metric used by BOSTIn.\n")
summ_list <- append(summ_list, paste("You can find a summary histogram of the nRCFV for each gene here:\n",summ_pdf,"\n"))
summ_list <- append(summ_list, "\n\nGenewise nRCFV summary statistics\n")

summ_list <- append(summ_list, paste("mean\t", summ_mean, "\n"))
summ_list <- append(summ_list, paste("median\t", summ_median, "\n"))
summ_list <- append(summ_list, paste("standard deviation\t", summ_sd, "\n"))
summ_list <- append(summ_list, paste("median absolute deviation\t", summ_mad, "\n"))
summ_list <- append(summ_list, paste("minimum\t", summ_min, "\n"))
summ_list <- append(summ_list, paste("maximum\t", summ_max, "\n"))
summ_list <- append(summ_list, "\n\nThe histogram is particularly useful for evaluating nRCFV. It shows the distribution of nRCFV values across the genes in your dataset, which, when compared to the average, can help you work out which taxa might end up causing problems in your final phylogenetic analysis. High nRCFV values indicate high levels of compositional heterogeneity.\nWe can also use the difference between the median and the mean of the tsnRCFV values to better understand what the distribution of compositional heterogeneity in your dataset looks like.\n")
summ_list <- append(summ_list, paste("The nRCFV values in your dataset range from ",summ_min," to ",summ_max,", with a mean of",summ_mean,"a standard deviation of",summ_sd," a median of ",summ_median," and a median absolute deviation of ",summ_mad,".\n"))
summ_list <- append(summ_list, "From the median, we can identify potentially problematic taxa using two metrics, the standard deviation and the median absolute deviation. In a normal distribution, 95% of data should exist within 2 standard deviations of the mean, but single instances of very large values can inflate the mean by quite a bit. If we instead measure 2 standard deviations from the median, we can more comfortably identify any potential outlier. The median absolute deviation, meanwhile, is a bit like the median's version of a standard deviation. It is the median of how far away each value in the dataset is from the median. It is a stricter measurement, and so all genes that are 2 standard deviations from the median are marked with a Red Flag so that you can seriously consider excluding them from your analysis, whereas all genes that are 2 median absolute deviations from the median are marked with a Yellow Flag to warn you that they might be a problem.\nYour Red Flag genes are:\n")
summ_list <- append(summ_list, paste(redList, sep="\n"))
summ_list <- append(summ_list, "\nYour Yellow Flag genes are:\n")
summ_list <- append(summ_list, paste(yellowList, sep="\n"))
cat(summ_list, file=summ_file)
