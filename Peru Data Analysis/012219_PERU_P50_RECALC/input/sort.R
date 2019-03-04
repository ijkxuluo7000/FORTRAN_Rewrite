setwd("C:/Users/wanju/Desktop/input")
getwd()

filelist = list.files(pattern = "txt")

out_file = vector()

for(i in filelist){
  
  the_file = readLines(i)
  the_file = c(i, the_file)
  the_file = c(the_file, "\n")
  out_file = c(out_file, the_file)
  
}

writeLines(out_file, con = "all_in_one.txt")