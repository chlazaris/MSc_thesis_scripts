#install the calibrate package
install.packages("calibrate")

#load calibrate package
library(calibrate)

setwd("/Users/chlazaris/Desktop/Desktop/MSC_lab/Results/Gal_BLAST_result_graphs")

Gal80_data <- read.csv("Sc_Gal80_BLAST_top_hits.csv", header = TRUE)
Gal80_data

x <- Gal80_data$Bit_score_1
y <- Gal80_data$Bit_score_2

# plot(Bit_score_1,Bit_score_2)
# xaxs and yaxs set to "i" to meet at zero
# limits for axes are set as well

#pdf("Sc_Gal80p_BLAST_top_hits.pdf")
plot(x,y, main="Gal80p top BLAST hit comparison", xlab="First_BLAST_hit_score", ylab="Second_BLAST_hit_score", xaxs="i", yaxs="i", xlim=c(0,1000), ylim=c(0,1000))
	 #use identify function to get species names that correspond
	 #to individual dots
	 identify(x,y,labels=Gal80_data$Species, offset=0.5, cex=0.7)
	 #add line of red color. a is for
	 #intercept and b i for slope
	 abline(a=0, b=1, col="red")
#dev.off()
