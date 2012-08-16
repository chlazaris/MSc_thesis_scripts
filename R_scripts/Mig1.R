#install the calibrate package
install.packages("calibrate")

#load calibrate package
library(calibrate)

setwd("/Users/chlazaris/Desktop/Desktop/MSC_lab/Results/Gal_BLAST_result_graphs")

Mig1_data <- read.csv("Sc_Mig1p_BLAST_top_hits.csv", header = TRUE)
Mig1_data

x <- Mig1_data$Bit_score_1
y <- Mig1_data$Bit_score_2

# plot(Bit_score_1,Bit_score_2)
# xaxs and yaxs set to "i" to meet at zero
# limits for axes are set as well

#pdf("Sc_Mig1p_BLAST_top_hits.pdf")
	plot(x,y, main="Mig1p top BLAST hit comparison", xlab="First_BLAST_hit_score", ylab="Second_BLAST_hit_score", xaxs="i", yaxs="i", xlim=c(0,1100), ylim=c(0,1100))
	identify(x,y,labels=Mig1_data$Species, offset=0.5, cex=0.7)
	#add line of red color. a is for
	#intercept and b i for slope
	abline(a=0, b=1, col="red")
#dev.off()
