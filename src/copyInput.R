### aSPIre ###
# description:  copy Skyline input to the .raw file location
# input:        sample list, Skyline inputs (.ssl, .fasta)
# output:       -
# author:       HPR

library(dplyr)


print("----------------------------------------")
print("2) COPY INPUT FILES TO RAW FILE LOCATION")
print("----------------------------------------")

protein_name = "IL37b"

### INPUT ###
sample_list = read.csv("data/sample_list.csv", stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name == protein_name, ]

### MAIN PART ###

# ----- check that all raw files are there -----
rawLoc = unique(sample_list$raw_file_path)
rawFiles = sapply(rawLoc, list.files, pattern = ".raw",full.names=T) %>%
  unlist()

if(all(sample_list$raw_file %in% basename(rawFiles))) {
  print("found all raw files :)")
  print("copying input to the raw file directory ...")
  
} else {
  
  rawLoc = gsub("/Volumes/DATA16040","/Volumes/FS/DATA16040",rawLoc)
  rawFiles = sapply(rawLoc, list.files, pattern = ".raw",full.names=T) %>%
    unlist()
  if(all(sample_list$raw_file %in% basename(rawFiles))) {
    print("found all raw files :)")
    print("copying input to the raw file directory ...")
  } else {
    print("sorry, the raw files are not in the specified directory ...")
  }
}

# ----- copy files -----

TOCOPY = c(paste0("results/",protein_name,"/",protein_name,".ssl"),
           paste0("results/",protein_name,"/",protein_name,".fasta"))

wd = getwd()
system(paste0("chmod -R 777 ", TOCOPY))
system(paste0("cd ",rawLoc[i],"; cp -rf ",wd,"/",TOCOPY," ./"))

# only copy if all .raw files are there
if(all(sample_list$raw_file %in% basename(rawFiles))) {
  
  able2copy = rep(F, length(TOCOPY))
  
  for (i in 1:length(rawLoc)) {
    able2copy[file.exists(paste0(rawLoc[i], basename(TOCOPY)))] = T
    
    while (! all(able2copy)) {
      able2copy = file.copy(from = TOCOPY[!able2copy], to = rawLoc[i],
                            copy.mode = F, copy.date = T)
    }
    
  }
  
  
}



file.copy("/Volumes/DATA16040/DATA/MS_data/MPI_BPC/ANALYSIS/WAITUCK/Mascot_distiller/aSyn/W_Soh_160821_170921_Goe_aSyn_A1_0h_R1.raw", "./")



ssl = read.table("results/IL37b/IL37b.ssl", stringsAsFactors = F, header = T, sep = "\t")
ssl$file = paste0("WAITUCK/Mascot_distiller/IL37b/",ssl$file)

write.table(ssl,"results/IL37b/IL37b.ssl",row.names = F,sep = "/t")


