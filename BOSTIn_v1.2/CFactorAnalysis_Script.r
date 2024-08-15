#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
#Take Command Line Arguments. Read the table in the supplied file and make the histogram of the table.
t_filename <- args[1]
tCdata <- read.table(t_filename, skip=1)
t_pdf <- paste0(tools::file_path_sans_ext(t_filename),".histogram.pdf")
pdf(file = t_pdf, width=4, height=4)
hist(tCdata$V4, main="Histogram of taxa C-Factor Score values", xlab="taxa C-Factor Score Values")
dev.off()

#Now do the useful statistics - First Mean, then median, then standard deviation, median absolute devitation, and percentage with a CFactor equal to or greater than 20, then 10.
t_summ_file <- paste0(tools::file_path_sans_ext(t_filename),".summary.txt")
t_summ_list <- c()
t_summ_list <- append(t_summ_list, "taxa C-Factor summary statistics\n")
t_mean = mean(tCdata$V4)
t_summ_list <- append(t_summ_list, paste("mean\t", t_mean, "\n"))
t_median = median(tCdata$V4)
t_summ_list <- append(t_summ_list, paste("median\t", t_median, "\n"))
t_sd = sd(tCdata$V4)
t_summ_list <- append(t_summ_list, paste("standard deviation\t", t_sd, "\n"))
t_mad = mad(tCdata$V4)
t_summ_list <- append(t_summ_list, paste("median absolute deviation\t", t_mad, "\n"))
cat(t_summ_list, file=t_summ_file)
t_cent20 = (length(tCdata$V4[tCdata$V4>=20])/length(tCdata$V4))*100
t_summ_list <- append(t_summ_list, paste("Percentage C-Factor >=20\t", t_cent20, "\n"))
cat(t_summ_list, file=t_summ_file)
t_cent10 = (length(tCdata$V4[tCdata$V4>=10])/length(tCdata$V4))*100
t_summ_list <- append(t_summ_list, paste("Percentage C-Factor >=10\t", t_cent10, "\n"))
cat(t_summ_list, file=t_summ_file)

#Print names of taxa, alongside their C Factors, with a C Factor of less than 10 to the red flags file.
redfile <- paste0(tools::file_path_sans_ext(t_filename),".redflags.txt")
t_RedFlags<-file(redfile)
redList <- c()
for(i in 1:length(tCdata$V4)){
	if(tCdata$V4[i] < 10){
		redTaxa <- paste(tCdata$V1[i], "\t", tCdata$V4[i], "\n")
		redList<- append(redList, redTaxa)
	}
}
cat(redList, file=t_RedFlags)

#Print names of taxa, alongside their C Factors, with a C Factor of less than 20 to the yellow flags file.
yellowfile <- paste0(tools::file_path_sans_ext(t_filename),".yellowflags.txt")
t_YellowFlags<-file(yellowfile)
yellowList <- c()
for(i in 1:length(tCdata$V4)){
	if(tCdata$V4[i] < 20){
		yellowTaxa <- paste(tCdata$V1[i], "\t", tCdata$V4[i], "\n")
		yellowList<- append(yellowList, yellowTaxa)
	}
}
cat(yellowList, file=t_YellowFlags)
