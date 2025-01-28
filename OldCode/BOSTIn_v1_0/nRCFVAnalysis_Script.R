#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
t_filename <- args[1]
tRCFVdata <- read.table(t_filename, skip=1)
t_pdf <- paste0(tools::file_path_sans_ext(t_filename),".histogram.pdf")
pdf(file = t_pdf, width=4, height=4)
hist(tRCFVdata$V2, main="Histogram of ntsRCFV values", xlab="ntsRCFV Values")
dev.off()

t_summ_file <- paste0(tools::file_path_sans_ext(t_filename),".summary.txt")
t_summ_list <- c()
t_summ_list <- append(t_summ_list, "ntsRCFV summary statistics\n")
t_mean = mean(tRCFVdata$V2)
t_summ_list <- append(t_summ_list, paste("mean\t", t_mean, "\n"))
t_median = median(tRCFVdata$V2)
t_summ_list <- append(t_summ_list, paste("median\t", t_median, "\n"))
t_sd = sd(tRCFVdata$V2)
t_summ_list <- append(t_summ_list, paste("standard deviation\t", t_sd, "\n"))
t_mad = mad(tRCFVdata$V2)
t_summ_list <- append(t_summ_list, paste("median absolute deviation\t", t_mad, "\n"))
cat(t_summ_list, file=t_summ_file)

yellowfile <- paste0(tools::file_path_sans_ext(t_filename),".CompHet.Taxa.yellowflags.txt")
t_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(tRCFVdata$V2)){
	if(tRCFVdata$V2[i] > (median(tRCFVdata$V2)+(2*mad(tRCFVdata$V2)))){
		yellowTaxa <- paste(tRCFVdata$V1[i], "\t", tRCFVdata$V2[i], "\n")
		yellowList<- append(yellowList, yellowTaxa)
	}
}
cat(yellowList, file=t_YellowFlags)

redfile <- paste0(tools::file_path_sans_ext(t_filename),".CompHet.Taxa.redflags.txt")
t_RedFlags<-file(redfile)
redList <- c()
for(i in 1:length(tRCFVdata$V2)){
	if(tRCFVdata$V2[i] > (median(tRCFVdata$V2)+(2*sd(tRCFVdata$V2)))){
		redTaxa <- paste(tRCFVdata$V1[i], "\t", tRCFVdata$V2[i], "\n")
		redList<- append(redList, redTaxa)
	}
}
cat(redList, file=t_RedFlags)

c_filename <- args[2]
cRCFVdata <- read.table(c_filename, skip=1)
c_pdf <- paste0(tools::file_path_sans_ext(c_filename),".histogram.pdf")
pdf(file = c_pdf, width=4, height=4)
hist(cRCFVdata$V2, main="Histogram of ncsRCFV values", xlab="ncsRCFV Values")
dev.off()

c_summ_file <- paste0(tools::file_path_sans_ext(c_filename),".summary.txt")
c_summ_list <- c()
c_summ_list <- append(c_summ_list, "ncsRCFV summary statistics\n")
c_mean = mean(cRCFVdata$V2)
c_summ_list <- append(c_summ_list, paste("mean\t", c_mean, "\n"))
c_median = median(cRCFVdata$V2)
c_summ_list <- append(c_summ_list, paste("median\t", c_median, "\n"))
c_sd = sd(cRCFVdata$V2)
c_summ_list <- append(c_summ_list, paste("standard deviation\t", c_sd, "\n"))
c_mad = mad(cRCFVdata$V2)
c_summ_list <- append(c_summ_list, paste("median absolute deviation\t", c_mad, "\n"))
cat(c_summ_list, file=c_summ_file)

c_yellowfile <- paste0(tools::file_path_sans_ext(c_filename),".CompHet.Character.yellowflags.txt")
c_YellowFlags<-file(c_yellowfile)
c_yellowList <- c()
for(i in 1:length(cRCFVdata$V2)){
	if(cRCFVdata$V2[i] > (median(cRCFVdata$V2)+(2*mad(cRCFVdata$V2)))){
		c_yellowTaxa <- paste(cRCFVdata$V1[i], "\t", cRCFVdata$V2[i], "\n")
		c_yellowList<- append(c_yellowList, c_yellowTaxa)
	}
}
cat(c_yellowList, file=c_YellowFlags)

redfile <- paste0(tools::file_path_sans_ext(c_filename),".CompHet.Character.redflags.txt")
c_RedFlags<-file(redfile)
redList <- c()
for(i in 1:length(cRCFVdata$V2)){
	if(cRCFVdata$V2[i] > (median(cRCFVdata$V2)+(2*sd(cRCFVdata$V2)))){
		redTaxa <- paste(cRCFVdata$V1[i], "\t", cRCFVdata$V2[i], "\n")
		redList<- append(redList, redTaxa)
	}
}
cat(redList, file=c_RedFlags)
