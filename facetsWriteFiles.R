#
# Take directory of FACETS objects and write to tsv
#

setwd(paste0("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")) 
source("helperFunctions.R")
source("facetsAnalysisLibrary.R")
library(stringr)
fileDf<- read.table("facetsRCancer/listR1File.txt", sep="\t")
fileList <- fileDf$V1
names(fileList) <- rownames(fileList)
reg <- "_([^/]+)/"

normal_ref <- 40
ref_match <- str_extract(fileList[[normal_ref]], reg)[1]
ref_sample <- substr(ref_match, 2, nchar(ref_match) - 1)

for(fileList.index in seq(normal_ref, length(fileList))){
  if(fileList.index == normal_ref) next
  target_match <- str_extract(fileList[[fileList.index]], reg)[1]
  target_sample <- substr(target_match, 2, nchar(target_match) - 1)
  target_class <- substr(target_sample, 2,2)
  
  
  res_dir <- if(target_class == "N") "facetsRNormal" else "facetsRCancer"
  xxFilename <- paste0(res_dir, "/facetsG5XX_", normal_ref, "_", fileList.index, ".rds")
  fitFilename <- paste0(res_dir, "/facetsG5Fit_", normal_ref, "_", fileList.index, ".rds")
  xx <- readRDS(xxFilename)
  fit <- readRDS(fitFilename)
  saveSnps <- xx$jointseg
  saveFit <- fit$cncf
  
  out_dir <- paste0("output/", "Sample_", target_sample)
  dir.create(file.path(out_dir), showWarnings = FALSE)
  xxOutputFilename <- paste0(out_dir, "/", target_sample, "--", ref_sample, ".procSample-jseg.cnv.facets.v0.5.2.txt")
  fitOutputFilename <- paste0(out_dir, "/", target_sample, "--", ref_sample, ".cnv.facets.v0.5.2.txt")
  write.table(saveSnps, file = xxOutputFilename, quote=TRUE, sep = "\t", row.names = FALSE  )
  write.table(saveFit, file = fitOutputFilename, quote=TRUE, sep = "\t", row.names = FALSE  )
}


