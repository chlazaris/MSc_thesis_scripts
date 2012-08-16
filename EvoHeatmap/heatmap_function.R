# R script that accepts .csv file (result of EvoHeatMap run)
# as input and returns a heatmap with clustering based on mean 
# row values.

# Type in terminal:
# RScript heatmap_function.R [name of .csv file]
# or use it as part of EvoHeatmap

# Used as source of help http://flowingdata.com/2010/01/21/how-to-make-a-heatmap-a-quick-and-easy-solution/

Args <- commandArgs(trailingOnly = TRUE);

data <- read.csv(Args[1])
name <- Args[1]

sorted_data <- data[order(data$Species),]

# define the row names
row.names(sorted_data) <- sorted_data$Species

# As you have set up the first column
# to be species names, we can get rid
# of the initial first column
ncol <- ncol(sorted_data)
sorted_data_2 <- sorted_data[,2:ncol]

# save as matrix in order to create
# the heatmap later on
g <- data.matrix(sorted_data_2)

#load required libraries
library("gplots")
library(dichromat)

#save as pdf with the selected dimensions
pdf(paste(name,".pdf",sep=""), height=8, width = 10)
heatmap.2(g, col=colorschemes$BrowntoBlue.10, scale="none", na.rm=TRUE,  
		  xlab="Proteins", ylab="Species", main=Args[2], key=TRUE, 
		  keysize=0.9, symkey=FALSE, #density.info="histogram",
		  #denscol="black", 
		  trace="none", cexCol=1.0, cexRow=0.6)
dev.off()