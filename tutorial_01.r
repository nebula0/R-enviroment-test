#========== Loading expression data
#==========

# show current directory on terminal
getwd(); 

# set workingDir to current directory
workingDir = "." 
setwd(workingDir)

# load package
library(WGCNA);

# don't known meaning
options(stringsAsFactors = FALSE);

# read data
femData = read.csv("LiverFemale3600.csv");

# dim similar to python's shape
dim(femData)
names(femData)

# remove auxiliary data transit data
# not quite understand
datExpr0 = as.data.frame(t(femData[, -c(1:8)])); names(datExpr0) = femData$substanceBXH;
rownames(datExpr0) = names(femData)[-c(1:8)];

#========== Checking data for excessive missing values 
#========== and identification of outlier microarray samples

gsg = goodSamplesGenes(datExpr0, verbose = 3); 
gsg$allO

if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed: if (sum(!gsg$goodGenes)>0) printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data: datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

R.version

Rhome 
