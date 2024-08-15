#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
t_filename <- args[1]
tDEdata <- read.table(t_filename, skip=1)
t_pdf <- paste0(tools::file_path_sans_ext(t_filename),".histogram.pdf")
pdf(file = t_pdf, width=4, height=4)
hist(tDEdata$V5, main="Histogram of t-DE Score values", xlab="t-DE Score Values")
dev.off()
#This section creates the summary file, then calculates the mean, median, standard deviation and median absolute deviation of the taxon specific DE Scores.
t_summ_file <- paste0(tools::file_path_sans_ext(t_filename),".summary.txt")
t_summ_list <- c()
t_summ_list <- append(t_summ_list, "t-DE summary statistics\n")
t_mean = mean(tDEdata$V5)
t_summ_list <- append(t_summ_list, paste("mean\t", t_mean, "\n"))
t_median = median(tDEdata$V5)
t_summ_list <- append(t_summ_list, paste("median\t", t_median, "\n"))
t_sd = sd(tDEdata$V5)
t_summ_list <- append(t_summ_list, paste("standard deviation\t", t_sd, "\n"))
t_mad = mad(tDEdata$V5)
t_summ_list <- append(t_summ_list, paste("median absolute deviation\t", t_mad, "\n"))
cat(t_summ_list, file=t_summ_file)

#Based on the MAD from the median, this identifies red flag taxa and places them in the output txt file.
redfile <- paste0(tools::file_path_sans_ext(t_filename),".SiteSat.Taxa.redflags.txt")
t_RedFlags<-file(redfile)
redList <- c()
for(i in 1:length(tDEdata$V5)){
	if(tDEdata$V5[i] < (median(tDEdata$V5)-(2*mad(tDEdata$V5)))){
		redTaxa <- paste(tDEdata$V1[i], "\t", tDEdata$V5[i], "\n")
		redList<- append(redList, redTaxa)
	}
}
cat(redList, file=t_RedFlags)

#Based on the standard deviation from the median, this identifies yellow flag taxa and places them in the output txt file.
yellowfile <- paste0(tools::file_path_sans_ext(t_filename),".SiteSat.Taxa.yellowflags.txt")
t_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(tDEdata$V5)){
	if(tDEdata$V5[i] < (median(tDEdata$V5)-(2*sd(tDEdata$V5)))){
		yellowTaxa <- paste(tDEdata$V1[i], "\t", tDEdata$V5[i], "\n")
		yellowList<- append(yellowList, yellowTaxa)
	}
}
cat(yellowList, file=t_YellowFlags)
