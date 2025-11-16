
################################################################################
#                                                                              #
#                R code for plotting the results of PG-SCUnK                   #
#                                                                              #
################################################################################

# USAGE:
# Rscript --vanilla PG-SCUnK_plot.R file.stats.txt out.pdf
# see PG-SCUnK webpage for details : https://github.com/cumtr/PG-SCUnK
# Developer : Tristan CUMER - t.cumer.sci[at]gmail.com

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
  install.packages("Ternary", repos = "http://cran.us.r-project.org", quiet = TRUE)
}

invisible(.libPaths())
library(Ternary)

#### read the data #####

Data = read.table(args[1], fill=TRUE)
Data = na.omit(Data)

#### make the plot #####

# ternary plots

pdf(file = paste0(args[2],".ternary.pdf"), width = 10, height = 10)

layout(matrix(1:4, ncol = 2, byrow = TRUE))
par(mar = rep(0.5, 4))

TernaryPlot(region = list(min=c(0,0,0), max=c(100,100,100)), alab = "Unique", blab = "Duplicated", clab = "Split")
legend("topleft", legend = "0-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

TernaryPlot(region = list(min=c(50,0,0), max=c(100,50,50)), alab = "Unique", blab = "Duplicated", clab = "Split")
legend("topleft", legend = "50-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

TernaryPlot(region = list(min=c(90,0,0), max=c(100,50,50)), alab = "Unique", blab = "Duplicated", clab = "Split")
legend("topleft", legend = "90-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

TernaryPlot(region = list(min=c(97.5,0,0), max=c(100,50,50)), alab = "Unique", blab = "Duplicated", clab = "Split")
legend("topleft", legend = "97.5-100%", bty = "n")
for(i in 1:nrow(Data)){
  AddToTernary(graphics::points, Data[i,c(2:4)], pch = 20)
}

dev.off()


# barplots

dat = Data[,c(3:5)]/rowSums(Data[,c(3:5)])
dat = as.matrix(dat)
rownames(dat) = Data[,1]
if(ncol(dat)>1){
    dat = t(dat); dat = dat[, rev(colnames(dat))]
}

pdf(file = paste0(args[2],"barplot.pdf"), width = 10, height = nrow(Data)+4)

par(mar = c(4,8,3,3))
barplot(dat, width = 1, space = 0, horiz = TRUE, las = 2)
text(0.02, (1:nrow(Data))-0.5, paste0("U: ", round(dat[1,], 3), " - D: ", round(dat[2,], 3), " - S: ", round(dat[3,], 3)), adj = c(0,0.5), col = "white")

dev.off()


#########################
