

#######################################
## FUNCTION TO READ FCS FILES 
######################################

read.fcs <- function(fcs.dir.path,
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
    for (file in included.fcs) {cat(file, "\n")}
    
    
  }else{warning("No files were processed")}

  cat("-------------------------", "\n")
  if (length(excluded.fcs)>0) {
    cat("The following files from the directory were excluded:", "\n")
    for (file in excluded.fcs) {cat(file, "\n")}
  }else{cat("all fcs files were included", "\n")}
  
  cat("-------------------------", "\n")
  if (length(missing.fcs)>0) {
    cat("The following files were not found in the fcs directory:", "\n")
    for (file in missing.fcs) {cat(file, "\n")}
  }else{cat("No missing files", "\n")}
  
  return(fs)
}




############################################
## FUNCTION TO FILTER AND PROCESS PARAMETERS 
############################################

process.param <- function(flowframe,
                          cofactor=1000,
                         label.sep="_"){
  
  # Obtain channels and labels from flowframe.
  param.df <- data.frame(channel = flowframe@parameters@data$name, 
                         label = as.character(flowframe@parameters@data$desc))
  
  #This will omit FSC, SSC, and Time parameters
  param.df <- na.omit(param.df[str_detect(param.df$channel, "^(?!FSC)(?!SSC).*-A$"),])

  # Detect the name of the marker in the label (i.e. excluding fluorochrome if it is also in the label)
  # If label separator is specified.
  if(!is.null(label.sep)){
    label.regex <- paste("(.*)[", label.sep, "](.*)", sep="")
  } else{
    label.regex <- "(.*)"
  }
  
  param.df$marker <- str_match(param.df$label, label.regex)[,2]
  param.df$marker[is.na(param.df$marker)] <- param.df$label[is.na(param.df$marker)]
  param.df$cofactor <- cofactor 
  
  # Show the available markers and their default cofactors.
  cat("copy, paste, and edit the cofactors as needed in the next section:", "\n")
  cat("---------------------------", "\n")
  for(marker in param.df$marker){
    cat(paste("'", marker, "'", "=", cofactor, ",", sep=""), "\n")
  }
  return(param.df)
}


