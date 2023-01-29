### aSPIre ###
# description:  filter final assignments and create Skyline input files
# input:        inSPIRE final assignment, inSPIRE features, sample list
# output:       .ssl file for Skyline input, formatted final assignments,
#               fasta file for Skyline
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
source("src/invitroSPI_utils.R")

protein_name = snakemake@params[["protein_name"]]

suppressWarnings(dir.create(paste0("results/",protein_name,"/")))
suppressWarnings(dir.create(paste0("results/",protein_name,"/plots/")))
sink(file = paste0("results/",protein_name,"/log.txt"), append = F, split = T)

print("----------------------------------------------")
print("PARSE INSPIRE RESULTS AND CREATE SKYLINE INPUT")
print("----------------------------------------------")

print(protein_name)

spAngleCutoff = snakemake@params[["spAngleCutoff"]]
qValCutoff = snakemake@params[["qValCutoff"]]
rtCutoff = snakemake@params[["rtCutoff"]]

paste0("spectral angle cut-off: ", spAngleCutoff) %>% print()
paste0("q-value cut-off: ", qValCutoff) %>% print()
paste0("retention time cut-off: ", rtCutoff) %>% print()

### INPUT ###
sample_list = read.csv(snakemake@input[["sample_list"]], stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name == protein_name, ]

finalAssignments = read.csv(paste0("data/inSPIRE/",sample_list$final_assignments[1]),
                            stringsAsFactors = F)


### MAIN PART ###
# ----- 1) parse and filter final assignments -----
# filtering
finalAssignments = finalAssignments[which(finalAssignments$qValue < qValCutoff & 
                                      finalAssignments$spectralAngle > spAngleCutoff &
                                      finalAssignments$deltaRT < rtCutoff), ]

# match with sample list to get info
sample_list = sample_list %>%
  select(protein_name, substrateID, substrateSeq, digestTime, biological_replicate, raw_file) %>%
  rename(source = raw_file)
sample_list$source = gsub(".raw","",sample_list$source)

finalAssignments = finalAssignments %>% rename(scanNum = scan, pepSeq = peptide, ionScore = engineScore)
DB = left_join(finalAssignments, sample_list)


# ----- 2) add substrate sequences -----
seqs = DB$substrateSeq %>% unique() %>% na.omit()
Dseqs = lapply(seqs, function(x){
  read.fasta(file=paste0("data/sequences/",x),seqtype = "AA", as.string=TRUE)[1] %>% unlist()
}) %>%
  plyr::ldply()
Dseqs$substrateSeq = seqs
names(Dseqs)[1] = "substrate"

# pre-processing
DB = left_join(DB, Dseqs)
DB$substrateSeq = DB$substrate
DB$substrate = NULL


# ----- 3) mapping -----
# DB$productType = str_extract_all(DB$proteins, "^[:alpha:]{3}", simplify = T)
DB$productType = "PCP"
DB$spliceType = NA
DB$positions = NA

ASSIGNMENTS = mapping(DB)
ASSIGNMENTS = apply(ASSIGNMENTS, 2, as.character) %>%
  as.data.frame()
ASSIGNMENTS$scanNum = as.numeric(ASSIGNMENTS$scanNum)
ASSIGNMENTS$charge = as.numeric(ASSIGNMENTS$charge)
if ("deltaRT" %in% names(ASSIGNMENTS)) {
  ASSIGNMENTS$deltaRT = as.numeric(ASSIGNMENTS$deltaRT)
}

k = which(str_detect(ASSIGNMENTS$positions, "^[:digit:]+_[:digit:]+_[:digit:]+_[:digit:]+"))
ASSIGNMENTS$productType[k] = "PSP"


# ----- 4) create .ssl table -----
# add charge and modifications
# SKYLINE = AllFeatures %>%
#   mutate(source = str_extract_all(AllFeatures$PSMId, "^[:graph:]+(?=_[:digit:]{3,}_[:alnum:]{5,}$)", simplify = T),
#          scanNum = str_extract_all(AllFeatures$PSMId, "[:digit:]{3,}(?=_[:alnum:]{5,}$)", simplify = T) %>% as.numeric()) %>%
#   filter(PSMId %in% paste0(ASSIGNMENTS$source,"_",ASSIGNMENTS$scanNum,"_",ASSIGNMENTS$modifiedSequence)) %>%
#   mutate(ID = paste0(source,"_",scanNum)) %>%
#   select(source,scanNum,ID,charge,deltaRT) %>%
#   right_join(ASSIGNMENTS) %>%
#   unique()

SKYLINE = ASSIGNMENTS %>%
  mutate(ID = paste0(source,"_",scanNum))

# create input
input = data.frame(file = paste0(SKYLINE$source,".raw"),
                   scan = SKYLINE$scanNum,
                   charge = SKYLINE$charge,
                   sequence = SKYLINE$modifiedSequence,  # !!!
                   score = SKYLINE$qValue,
                   modifications = SKYLINE$modifiedSequence)

# ----- 5) create .fasta file -----
# pass only unique peptides
P = SKYLINE %>%
  distinct(pepSeq, .keep_all = T)
prots = P$pepSeq
protList = lapply(prots, function(x){x})
nm = paste0(P$source,"_",P$scanNum)
names(protList) = nm

# ----- 6) inform about retention time -----
print("retention time deviation:")
summary(SKYLINE$deltaRT/60) %>%
  print()


### OUTPUT ###
write.csv(SKYLINE, file = unlist(snakemake@output[["assignments_csv"]]),
          row.names = F)
save(SKYLINE, file = unlist(snakemake@output[["assignments_bin"]]))

write.table(input, file = unlist(snakemake@output[["ssl"]]),
            sep = "\t", row.names = F)
write.fasta(sequences = protList, names =nm,
            file.out = unlist(snakemake@output[["fasta"]]),as.string = T)


