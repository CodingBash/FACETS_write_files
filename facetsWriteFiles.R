#
# Take directory of FACETS objects and write to tsv
#

setwd(paste0("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")) 
source("helperFunctions.R")
source("facetsAnalysisLibrary.R")
library(stringr)

xxNew <- readRDS("./resources/xxJoinsegNew.rds")

fileDf<- read.table("./resources/listR1File.txt", sep="\t")
fileList <- fileDf$V1
names(fileList) <- rownames(fileList)

#
# Determine reference sample
#
reg <- "_([^/]+)/"
normal_ref <- 41 # This is the reference sample chosen
ref_match <- str_extract(fileList[[normal_ref]], reg)[1]
ref_sample <- substr(ref_match, 2, nchar(ref_match) - 1)

for(fileList.index in seq(1, length(fileList))){
  #
  # Pre-work
  #
  if(fileList.index == normal_ref) next
  target_match <- str_extract(fileList[[fileList.index]], reg)[1]
  target_sample <- substr(target_match, 2, nchar(target_match) - 1)
  target_class <- substr(target_sample, 2,2)
  
  res_dir <- if(target_class == "N") "facetsRNormal" else "facetsRCancer"
  res_dir <- paste0("resources/", res_dir)

  #
  # Load SNPs
  #
  saveSnps <- NA
  if(target_class == "N"){
    xxFilename <- paste0(res_dir, "/facetsG5XX_", normal_ref, "_", fileList.index, ".rds")
    xx <- readRDS(xxFilename)  
    saveSnps <- xx$jointseg
  } else if (target_class %in% c("T", "F", "M")){
    saveSnps <- xxNew[[paste0("n", normal_ref, "_", fileList.index)]]    
  } else {
    print(paste0("WARNING: Target class not recognized: ", target_class))
    return()
  }
  
  #
  # Load segment fit
  #
  fitFilename <- paste0(res_dir, "/facetsG5Fit_", normal_ref, "_", fileList.index, ".rds")
  fit <- readRDS(fitFilename)
  saveFit <- fit$cncf
  
  #
  # Save to directory
  #
  out_dir <- paste0("output/", "Sample_", target_sample)
  dir.create(file.path(out_dir), showWarnings = FALSE)
  xxOutputFilename <- paste0(out_dir, "/", target_sample, "--", ref_sample, ".procSample-jseg.cnv.facets.v0.5.2.txt")
  fitOutputFilename <- paste0(out_dir, "/", target_sample, "--", ref_sample, ".cnv.facets.v0.5.2.txt")
  write.table(saveSnps, file = xxOutputFilename, quote=TRUE, sep = "\t", row.names = FALSE  )
  write.table(saveFit, file = fitOutputFilename, quote=TRUE, sep = "\t", row.names = FALSE  )
}

