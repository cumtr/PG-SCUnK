
################################################################################
#                                                                              #
#                R code for plotting the results of PG-SCUnK                   #
#                                                                              #
#                 Tristan CUMER - t.cumer.sci[at]gmail.com                     # 
#                                                                              #
################################################################################

# USAGE : 
# Rscript --vanilla PG-SCUnK_plot.R file.stats.txt out.pdf

#### Set up the environement #####

# load the arguements provided
args = commandArgs(trailingOnly=TRUE)

# check if the arguements match the expectations
if (length(args)==0) {
  stop("Usage : Rscript --vanilla PG-SCUnK_plot.R file.stats.txt out.pdf \n At least one argument must be supplied (output from PG-SCUnK, file.stats.txt).")
} else if (length(args)==1) {
  # default output file
  args[2] = "out.pdf"
}

# Check if package missing, if yes install, then load it
if(!require(Ternary)){
  install.packages("Ternary",repos = "http://cran.us.r-project.org")
}
library(Ternary)

#### read the data #####

Data = read.table(args[1], fill=TRUE)
Data = na.omit(Data)

#### make the plot #####

pdf(file = args[2], width = 10, height = 10)

layout(matrix(1:4, ncol = 2, byrow = TRUE))
par(mar = rep(0.5, 4))

TernaryPlot(region = list(min=c(0,0,0), max=c(100,100,100)), alab = "Unique", blab = "Duplicated", clab = "Colapsed")
legend("topleft", legend = "0-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

TernaryPlot(region = list(min=c(50,0,0), max=c(100,50,50)), alab = "Unique", blab = "Duplicated", clab = "Colapsed")
legend("topleft", legend = "50-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

TernaryPlot(region = list(min=c(90,0,0), max=c(100,50,50)), alab = "Unique", blab = "Duplicated", clab = "Colapsed")
legend("topleft", legend = "90-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

TernaryPlot(region = list(min=c(97.5,0,0), max=c(100,50,50)), alab = "Unique", blab = "Duplicated", clab = "Colapsed")
legend("topleft", legend = "97.5-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

dev.off()

#########################
