

#######################################
## FUNCTION TO READ FCS FILES 
######################################

read.fsc <- function(fcs.dir.path,
                     metadata,
                     transform=FALSE,
                     downsample=NULL,
                     truncate.max=FALSE) {
  
  metadata.files <- paste(metadata[,1], ".fcs", sep="") 
  all.fcs <- list.files(path=fcs.dir.path, pattern = "\\.fcs")
  
  included.fcs <- metadata.files[metadata.files %in% all.fcs]
  excluded.fcs <- all.fcs[!(all.fcs %in% metadata.files)]
  missing.fcs <- metadata.files[!(metadata.files %in% all.fcs)]
  
  if (length(included.fcs)>0) {
   fs <- read.flowSet(files=included.fcs, 
                      path = fcs.dir.path,
                      transformation=transform,
                      which.lines=downsample,
                      truncate_max_range = truncate.max)
  
    cat("The following files were included:", "\n")
    for (file in included.fcs) {cat(file)}
  }else{warning("No files were processed")}
  
  if (length(excluded.fcs)>0) {
    cat("The following files from the directory were excluded:", "\n")
    for (file in excluded.fcs) {cat(file)}
  }else{cat("all fcs files were included", "\n")}
  
  if (length(missing.fcs)>0) {
    cat("The following files were not found in the fcs directory:", "\n")
    for (file in missing.fcs) {cat(file)}
  }else{cat("No missing files", "\n")}
  
}