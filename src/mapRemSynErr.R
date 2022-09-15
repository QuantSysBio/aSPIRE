### aSPIre ###
# description:  map quantified peptides to assign splice types and remove synthesis errors
# input:        raw intensities, aggregated over biological and technical replicates
# output:       filtered raw intensities
# author:       HPR

library(dplyr, warn.conflicts = F)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/aSPIre_plotting.R")

protein_name = snakemake@params[["protein_name"]]
sink(file = paste0("results/",protein_name,"/log.txt"), append = T, split = T)

print("----------------------------------------")
print("MAP PEPTIDES AND REMOVE SYNTHESIS ERRORS")
print("----------------------------------------")


### INPUT ###
load(snakemake@input[["quantities_raw"]])
load(snakemake@input[["assignments_bin"]])

sample_list = read.csv(snakemake@input[["sample_list"]], stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name == protein_name, ]


### MAIN PART ###
# ----- map and assign splice types -----
QUANTITIES = as.data.frame(QUANTITIES)
QUANTITIES$productType = "PCP"
QUANTITIES$spliceType = NA
QUANTITIES$positions = NA

DB = mapping(QUANTITIES)
DB$productType[!is.na(DB$spliceType)] = "PSP"


# ----- retrieve peptides that were assigned at t=0 -----
# remove any whitespaces
SKYLINE$digestTime = str_extract_all(SKYLINE$digestTime, pattern = "[:digit:]+", simplify = T) %>% as.numeric()
DB$digestTime = str_extract_all(DB$digestTime, pattern = "[:digit:]+", simplify = T) %>% as.numeric()

SynErrors = SKYLINE[SKYLINE$digestTime %in% c(0,"CTRL","0"),]
errors = gsub("I","L",SynErrors$pepSeq) %>% unique()

paste0("peptides that were identified at t=0: ", length(errors)) %>%
  print()

# PSP synthesis errors based on sequence
if (length(errors) > 0) {
  
  ErrPeps = c()
  print("removing PSP synthesis errors based on MS2 identification...")
  
  # remove exact matches to synthesis errors
  pepU = gsub("I","L",DB$pepSeq[DB$productType == "PSP"]) %>% unique()
  k = which(pepU %in% errors)
  if (length(k) > 0) {
    pepRem = pepU[k]
    paste0("exact matches: ", length(pepRem)) %>% print()
    ErrPeps = c(ErrPeps,pepRem)
  }
  
  # remove precursors
  psp_remove = sapply(pepU, function(x){
    any(grepl(x,errors))
  }) %>%
    which()
  
  if (length(psp_remove) > 0) {
    pepRem = pepU[psp_remove]
    paste0("precursors: ", length(pepRem)) %>% print()
    ErrPeps = c(ErrPeps,pepRem)
  }
  
  
  # PCP synthesis errors based on kinetic
  print("removing PCP synthesis errors based on MS1 intensity...")
  pepC = gsub("I","L",DB$pepSeq[DB$productType == "PCP"]) %>% unique()
  j = which(pepC %in% errors)
  if (length(j) > 0) {
    pepQM = pepC[j]
    paste0(length(pepQM)," exact matches are being evaluated for their kinetic") %>% print()
    
    QMtable = DB %>%
      filter(gsub("I","L",pepSeq) %in% pepQM) %>%
      group_by(pepSeq,biological_replicate,digestTime) %>%
      arrange(digestTime, .by_group = T) %>%
      mutate(technical_replicate = row_number()) %>%
      ungroup() %>%
      group_by(pepSeq,biological_replicate,technical_replicate) %>%
      summarise(int_pasted = paste(intensity, collapse = ";"),
                times_pasted = paste(digestTime, collapse = ";"))
    
    tps = strsplit(QMtable$times_pasted[1],";") %>% unlist() %>% as.numeric() %>% sort()
    intens = str_split_fixed(QMtable$int_pasted,coll(";"),Inf)
    intens = apply(intens,2,as.numeric)
    
    pcp_rem = apply(intens,1,function(x){
      any(x[tps==0] > x[tps>0])
    }) %>%
      which()
   
    if(length(pcp_rem) > 0) {
      pepRem = gsub("I","L",QMtable$pepSeq[pcp_rem]) %>% unique()
      paste0("removing PCPs based on kinetic: ",length(pepRem)) %>%
        print()
      ErrPeps = c(ErrPeps,pepRem)
    }
    
  }
  
  # plot synthesis errors
  if (length(ErrPeps) > 0) {
    plotKinetics(DB[gsub("I","L",DB$pepSeq) %in% ErrPeps,],
                 outfile = paste0("results/",protein_name,"/plots/rawIntensities_SynErr.pdf"),
                 meanTech = F, earlyOnly = F, sortByInt = T)
  }
  
  # filter DB
  DB = DB[!gsub("I","L",DB$pepSeq) %in% gsub("I","L",ErrPeps),]
  
} else {
  print("no synthesis errors are being removed")
}

# ----- format intensity tables -----

QUANTfiltered = DB %>%
  group_by(pepSeq,biological_replicate,digestTime) %>%
  arrange(digestTime, .by_group = T) %>%
  mutate(technical_replicate = row_number()) %>%
  ungroup() %>%
  group_by(substrateID,pepSeq,productType,spliceType,positions,biological_replicate,technical_replicate) %>%
  summarise(int_pasted = paste(intensity, collapse = ";"),
            times_pasted = paste(digestTime, collapse = ";"))


# plot raw intensities
plotKinetics(QUANTfiltered,
             outfile = unlist(snakemake@output[["plot_raw"]]),
             meanTech = F, earlyOnly = F, sortByInt = T)

### OUTPUT ###
save(QUANTfiltered, file = unlist(snakemake@output[["quantities_filtered"]]))


