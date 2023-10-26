library(data.table)
library(vroom)

args = commandArgs(trailingOnly=TRUE)

path = args[1]
output_name = args[2]

files = list.files(path, full.names = TRUE) 
output = as.data.table(vroom(files))

fwrite(output, output_name, sep = '\t')
