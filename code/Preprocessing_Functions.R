

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
    
    return(fs)
    
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




#######################################
## FUNCTION TO FILTER PARAMETERS 
######################################

filter.param <- function(param.df,
                         AutoF.ch=NULL,
                         lab.sep="_",
                         cofactor=1000){
  
  #This will omit FSC, SSC, Time, and autofluorescence parameters
  param.df <- na.omit(param.df[str_detect(param.df$channel, "^(?!FSC)(?!SSC).*-A$"),])
  
  if(!(is.null(AutoF.ch))&(AutoF.ch %in% param.df$channel)){
    autof.regex <- paste("^(?!\\", AutoF.ch, ")", sep="")
    param.df <- na.omit(param.df[str_detect(param.df$channel, autof.regex),])
  }
  
  if(!is.null(lab.sep)){
    label.regex <- paste("(.*)[", lab.sep, "](.*)", sep="")
  } else{
    label.regex <- "(.*)"
  }
  
  # Detect the name of the marker in the label (i.e. excluding fluorochrome if it is also in the label)
  param.df$marker <- str_match(param.df$label, label.regex)[,2] 
  param.df$cofactor <- cofactor 
  
  return(param.df)
  for(marker in param.df$marker){
    cat(paste("'", marker, "'", "=", cofactor, ",", sep=""), "\n")
  }
}


