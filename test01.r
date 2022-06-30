getwd();
workingDir = ".";
library(WGCNA);
options(stringsAsFactors = FALSE);

femData = read.csv("LiverFemale3600.csv");
dim(femData);
names(femData)

