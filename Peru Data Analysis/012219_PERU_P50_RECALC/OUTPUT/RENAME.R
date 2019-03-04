# setwd("C:\\dos_exe\\OUTPUT")
# 
# filelist = list.files(pattern = "TXT")[1]
# 
# txt = readLines(filelist)
# whole_txt = vector()
# for(i in 1:length(txt)){
#   whole_txt = paste0(whole_txt, txt[i])
# }
# temp = unlist(strsplit(whole_txt, "AT A P50 OF  "))

rm(list = ls())

library(knitr)

get_directory = function(){
  args <- commandArgs(trailingOnly = FALSE)
  file <- "--file="
  rstudio <- "RStudio"
  
  match <- grep(rstudio, args)
  if(length(match) > 0){
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }else{
    match <- grep(file, args)
    if (length(match) > 0) {
      return(dirname(normalizePath(sub(file, "", args[match]))))
    }else{
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
  }
}

wd = get_directory()
setwd(wd)

get_p50 = function(FILEDIR){
  text <- readLines(FILEDIR,encoding="UTF-8")[25:250];
  num = as.numeric(unlist(lapply(text,FUN = strsplit,split = "          ")));
  return(num[which(num == min(num[seq(2,length(num),2)]))-1])
}

filelist = list.files(pattern = "TXT")

p50_list = vector()

for(i in 1:length(filelist)){
  
  p50_list[i] = get_p50(FILEDIR = filelist[i])
  
}

cat("\n")
cat("P50 from the normal approach")
cat("\n")
print(kable(as.data.frame(cbind(filelist, p50_list))))

get_p50_hill = function(FILEDIR){
  return(as.numeric(unlist(strsplit(readLines(FILEDIR)[length(readLines(FILEDIR))], "P50 = "))[2]))
}

p50_hill_list = vector()
for(i in 1:length(filelist)){
  
  p50_hill_list[i] = get_p50_hill(FILEDIR = filelist[i])
  
}

cat("\n")
cat("P50 from the hill approach")
cat("\n")
print(kable(as.data.frame(cbind(filelist, p50_hill_list))))
