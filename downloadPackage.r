getwd();
workingDir = ".";

orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("GO.db", "preprocessCore", "impute"));