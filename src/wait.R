### aSPIre ###
# description:  wait for Skyline results
# input:        -
# output:       -
# author:       HPR

library(dplyr)
protein_name = snakemake@params[["protein_name"]]
sink(file = paste0("results/",protein_name,"/log.txt"), append = T, split = T)

print("-----------------------")
print("WAIT FOR SKYLINE OUTPUT")
print("-----------------------")

### INPUT ###
ssl = snakemake@input[["ssl"]]
fasta = snakemake@input[["fasta"]]
result_file = unlist(snakemake@output[["skyline"]])

### MAIN PART ###
# ----- user info -----
print("------------------------------------------------------------------")
print("please copy the following files to your .raw file location:")
print(ssl)
print(fasta)
paste0("run Skyline and copy the following report to: results/",protein_name,": ")
print(basename(result_file))
print("------------------------------------------------------------------")

# ----- wait -----
while(!file.exists(result_file)) {
  
  paste0(Sys.time()," - waiting for Skyline output..") %>% print()
  Sys.sleep(30)
  
}

