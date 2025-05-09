#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
SAT <- read.table(filename, sep="\t", header=TRUE)

#Make the histogram
summ_pdf <- paste0(tools::file_path_sans_ext(filename),".histogram.pdf")
pdf(file = summ_pdf, width=4, height=4)
hist(SAT[,2], main="Histogram of t-DE Score values", xlab="t-DE Score Values")
dev.off()

#Make the Red and Yellow Lists
redfile <- paste0(tools::file_path_sans_ext(filename),".SiteSat.GeneWise.redflags.txt")
summ_RedFlags<-file(redfile)
redList <- c()
yellowfile <- paste0(tools::file_path_sans_ext(filename),".SiteSat.GeneWise.yellowflags.txt")
summ_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(SAT[,2])){
	if(SAT[i,2] <= 0){
		redGene <- paste(SAT[i,1], "\t", SAT[i,2], "\n")
		redList<- append(redList, redGene)
		yellowGene <- paste(SAT[i,1], "\t", SAT[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}
	else if(SAT[i,2] <= 0.354){
		yellowGene <- paste(SAT[i,1], "\t", SAT[i,2], "\n")
		yellowList<- append(yellowList, yellowGene)
	}

}
cat(redList, file=summ_RedFlags)
cat(yellowList, file=summ_YellowFlags)
#Write the Narrative Summary
summ_file <- paste0(tools::file_path_sans_ext(filename),".NarrativeSummary.txt")
summ_list <- c()

summ_list <- append(summ_list, "Welcome to the GeneWise batch selection narrative summary for the DE-Score, the amino acid saturation metric used by BOSTIn.\n")
summ_list <- append(summ_list, paste("You can find a summary histogram of each gene's whole dataset DE-Score here:\n",summ_pdf,"\n"))
summ_list <- append(summ_list, "\n\nGenewise DE summary statistics\n")
summ_mean = mean(SAT[,2])
summ_list <- append(summ_list, paste("mean\t", summ_mean, "\n"))
summ_median = median(SAT[,2])
summ_list <- append(summ_list, paste("median\t", summ_median, "\n"))
summ_sd = sd(SAT[,2])
summ_list <- append(summ_list, paste("standard deviation\t", summ_sd, "\n"))
summ_mad = mad(SAT[,2])
summ_list <- append(summ_list, paste("median absolute deviation\t", summ_mad, "\n"))
summ_min = min(SAT[,2])
summ_list <- append(summ_list, paste("minimum\t", summ_min, "\n"))
summ_max = max(SAT[,2])
summ_list <- append(summ_list, paste("maximum\t", summ_max, "\n"))
summ_list <- append(summ_list, "\n\nFor the DE-Score, values below 0 are considered to be saturated to an uninformative degree. We mark these genes as Red Flags. Since the uninformative ratio of the DE-Score is 0.177, we use two steps of uninformativeness and mark all genes with a DE-Score below 0.354 as Yellow Flags.\n Red Flags are sequences that may be too saturated to be informative, while Yellow Flags are ones to keep an eye on, as they may potentially cause problems in the dataset.\n")
summ_list <- append(summ_list, paste("You can find your list of Yellow Flag and Red Flag genes here:\n",yellowfile,"\n",redfile,"\nIn addition, we have included them below for convenience.\nYour Red Flag genes are:\n"))
summ_list <- append(summ_list, paste(redList, sep="\n"))
summ_list <- append(summ_list, "\nYour Yellow Flag genes are:\n")
summ_list <- append(summ_list, paste(yellowList, sep="\n"))
cat(summ_list, file=summ_file)
