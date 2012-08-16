#install the calibrate package
install.packages("calibrate")

#load calibrate package
library(calibrate)

setwd("/Users/chlazaris/Desktop/Desktop/MSC_lab/Results/Gal_BLAST_result_graphs")

Gal10_data <- read.csv("Sc_Gal10_BLAST_top_hits.csv", header = TRUE)
Gal10_data

x <- Gal10_data$Bit_score_1
y <- Gal10_data$Bit_score_2

# plot(Bit_score_1,Bit_score_2)
# xaxs and yaxs set to "i" to meet at zero
# limits for axes are set as well

#pdf("Sc_Gal10p_BLAST_top_hits.pdf")
	plot(x,y, main="Gal10p top BLAST hit comparison", xlab="First_BLAST_hit_score", ylab="Second_BLAST_hit_score", xaxs="i", yaxs="i", xlim=c(0,1500), ylim=c(0,1500))
	identify(x,y,labels=Gal10_data$Species, offset=0.5, cex=0.7)
	#add line of red color. a is for
	#intercept and b i for slope
	abline(a=0, b=1, col="red")
#dev.off()
