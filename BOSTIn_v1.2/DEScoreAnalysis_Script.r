#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
t_filename <- args[1]
tDEdata <- read.table(t_filename, skip=1)
t_pdf <- paste0(tools::file_path_sans_ext(t_filename),".histogram.pdf")
pdf(file = t_pdf, width=4, height=4)
hist(tDEdata$V4, main="Histogram of t-DE Score values", xlab="t-DE Score Values")
dev.off()

t_summ_file <- paste0(tools::file_path_sans_ext(t_filename),".summary.txt")
t_summ_list <- c()
t_summ_list <- append(t_summ_list, "t-DE summary statistics\n")
t_mean = mean(tDEdata$V4)
t_summ_list <- append(t_summ_list, paste("mean\t", t_mean, "\n"))
t_median = median(tDEdata$V4)
t_summ_list <- append(t_summ_list, paste("median\t", t_median, "\n"))
t_sd = sd(tDEdata$V4)
t_summ_list <- append(t_summ_list, paste("standard deviation\t", t_sd, "\n"))
t_mad = mad(tDEdata$V4)
t_summ_list <- append(t_summ_list, paste("median absolute deviation\t", t_mad, "\n"))
cat(t_summ_list, file=t_summ_file)

redfile <- paste0(tools::file_path_sans_ext(t_filename),".SiteSat.Taxa.redflags.txt")
t_RedFlags<-file(redfile)
redList <- c()
for(i in 1:length(tDEdata$V4)){
	if(tDEdata$V4[i] < 0){
		redTaxa <- paste(tDEdata$V1[i], "\t", tDEdata$V4[i], "\n")
		redList<- append(redList, redTaxa)
	}
}
cat(redList, file=t_RedFlags)

yellowfile <- paste0(tools::file_path_sans_ext(t_filename),".SiteSat.Taxa.yellowflags.txt")
t_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(tDEdata$V4)){
	if(tDEdata$V4[i] < 0.354){
		yellowTaxa <- paste(tDEdata$V1[i], "\t", tDEdata$V4[i], "\n")
		yellowList<- append(yellowList, yellowTaxa)
	}
}
cat(yellowList, file=t_YellowFlags)
