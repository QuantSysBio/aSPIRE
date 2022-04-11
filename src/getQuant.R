### aSPIre ###
# description:  read in Skyline output and parse relevant info
# input:        Skyline output (MS1_HPR.csv)
# output:       quantification for all peptides/time points
# author:       HPR

library(dplyr)
library(stringr)

protein_name = snakemake@params[["protein_name"]]
sink(file = paste0("results/",protein_name,"/log.txt"), append = T, split = T)

print("----------------------------------------")
print("RETRIEVE INTENSITIES FOR UNIQUE PEPTIDES")
print("----------------------------------------")


### INPUT ###
quant = read.csv(snakemake@input[["skyline"]], stringsAsFactors = F)
load(snakemake@input[["assignments_bin"]])
sample_list = read.csv(snakemake@input[["sample_list"]], stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name == protein_name, ]


### MAIN PART ###
# ----- remove all peptides that were not quantified -----
quantAssign = SKYLINE %>%
  select(-accessionGroup, -proteins) %>%
  mutate(digestTime = as.numeric(digestTime))

# print some info
NuniquePep = SKYLINE$pepSeq %>% unique() %>% length()
NquantPep = quant$Peptide.Sequence %>% unique() %>% length()

paste0(NquantPep," out of ",NuniquePep," unique peptides could be quantified") %>%
  print()


not = unique(SKYLINE$modifiedSequence[!SKYLINE$pepSeq %in% c(quant$Peptide.Sequence %>% unique())])
write.csv(SKYLINE[SKYLINE$modifiedSequence %in% not, ] %>% distinct(ID,.keep_all = T),
          paste0("results/",protein_name,"/notQuant.csv"),row.names = F)

# sanity check
# y = quant[quant$Peptide.Sequence == "AAHPGWF",]
# table(y$Product.Charge, y$Replicate.Name)

# ----- re-format quantification results -----
QUANT = quant %>%
  rename(ID = Protein.Name,
         source = File.Name,
         pepSeq = Peptide.Sequence,
         intensity = Normalized.Area,
         fragment = Fragment.Ion,
         charge = Product.Charge,
         mz = Product.Mz,
         RT = Average.Measured.Retention.Time,
         PTMs = Peptide.Modified.Sequence) %>%
  select(ID,source,pepSeq,intensity,fragment,charge,mz,RT,PTMs) %>%
  mutate(source = gsub(".raw","",source))

# replace NA intensities with 0
QUANT$intensity[QUANT$intensity == "#N/A"] = 0

# join with inSPIRE assignments
JOINED = left_join(QUANT, quantAssign, by = c("ID")) %>%
  filter(!is.na(spectralAngle))

JOINED = JOINED %>%
  rename(source = source.x,
         pepSeq = pepSeq.x,
         charge = charge.x) %>%
  select(-source.y, -pepSeq.y, -charge.y) %>%
  select(-digestTime,-biological_replicate)

# add info for each raw file using sample list
X = sample_list %>%
  mutate(source = gsub(".raw","",raw_file)) %>%
  select(source, digestTime, biological_replicate) %>%
  right_join(JOINED)

# ----- sum intensities for each peptides
X$intensity = as.numeric(X$intensity)

# all peptides have one precursor peak
# summed over all charges
Y = X %>%
  filter(fragment == "precursor") %>%
  group_by(pepSeq,source) %>% 
  mutate(intensity = sum(intensity)) %>%
  select(source,substrateID,digestTime,biological_replicate,pepSeq,intensity,substrateSeq,ID) %>%
  unique()

# Y %>%
#   group_by(pepSeq) %>%
#   summarise(n = length(unique(ID))) %>% View()


# summarise IDs
Z = Y %>%
  group_by(source,substrateID,digestTime,biological_replicate,pepSeq,intensity,substrateSeq) %>%
  summarise(assignments = paste(ID,collapse = ";"))

# ----- aggregate kinetics for replicates -----
QUANTITIES = Z %>%
  group_by(pepSeq,digestTime,biological_replicate) %>%
  mutate(mean_techRep = if (all(intensity == 0) & digestTime != 0) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T)) %>%
  ungroup() %>%
  group_by(pepSeq,digestTime) %>%
  mutate(mean_bioRep = if (all(intensity == 0) & digestTime != 0) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T))


### OUTPUT ###
save(QUANTITIES, file = unlist(snakemake@output[["quantities_raw"]]))

