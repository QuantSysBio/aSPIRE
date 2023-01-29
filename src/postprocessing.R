### aSPIre ###
# description:  normalising of intensities, plotting
# input:        filtered raw intensities, inSPIRE assignments
# output:       finalKinetics, data frame with normalised intensities
# author:       HPR

library(dplyr, warn.conflicts = F)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/aSPIre_plotting.R")

protein_name = snakemake@params[["protein_name"]]
sink(file = paste0("results/",protein_name,"/log.txt"), append = T, split = T)

print("--------------------------")
print("POSTPROCESSING OF KINETICS")
print("--------------------------")

### INPUT ###
load(snakemake@input[["quantities_filtered"]])
load(snakemake@input[["assignments_bin"]])


### MAIN PART ###
# ----- filter and plot final kinetics -----
print("normalisation...")

Q = QUANTfiltered %>%
  tidyr::separate_rows(int_pasted, times_pasted, sep = ";") %>%
  rename(intensity = int_pasted,
         digestTime = times_pasted) %>%
  mutate(digestTime = as.numeric(digestTime),
         intensity = as.numeric(intensity)) %>%
  as.data.frame()

# if there are no 0 hrs assignments
if (all(Q$digestTime > 0)) {
  Qmock = Q %>%
    ungroup() %>%
    distinct(pepSeq, biological_replicate, technical_replicate, .keep_all = T) %>%
    mutate(digestTime = 0,
           intensity = 0)
  
  Q = rbind(Q, Qmock) %>%
    as.data.frame() %>%
    arrange(digestTime, pepSeq)
  
} else if (length(unique(Q$biological_replicate[Q$digestTime == 0])) < length(unique(Q$biological_replicate))) {
 
  missing0 = Q %>%
    group_by(biological_replicate) %>%
    summarise(missing0 = all(digestTime > 0))
  missing0 = missing0$biological_replicate[missing0$missing0]
  
  Qmock = Q %>%
    ungroup() %>%
    filter(biological_replicate %in% missing0) %>%
    distinct(pepSeq, biological_replicate, technical_replicate, .keep_all = T) %>%
    mutate(digestTime = 0,
           intensity = 0)
  
  Q = rbind(Q, Qmock) %>%
    as.data.frame() %>%
    arrange(digestTime, pepSeq)
  
}


# --- aggregate intensities from I/L redundant peptides
Q = Q %>%
  mutate(pepSeq_IL = gsub("I","L",pepSeq)) %>%
  group_by(pepSeq_IL, biological_replicate, technical_replicate, digestTime) %>%
  summarise(intensity = sum(intensity)) %>%
  left_join(Q %>% mutate(pepSeq_IL = gsub("I","L",pepSeq)) %>% select(-intensity)) %>%
  ungroup() %>%
  select(-pepSeq_IL) %>%
  mutate(technical_replicate = ifelse(technical_replicate == 3, 1, technical_replicate),  # hard-coded for 2 technical replicates
         technical_replicate = ifelse(technical_replicate == 4, 2, technical_replicate))

# --- outlier filtering
# determine background from 0 hrs
t0 = Q %>%
  filter(digestTime == 0) %>%
  group_by(pepSeq, biological_replicate, technical_replicate) %>%
  summarise(max0 = max(intensity)) %>%
  ungroup()

# determine global noise level in digestion from distribution of 0 hrs values
noiseLevel = quantile(t0$max0, 0.5)

# set all values NA that are smaller than the background or the global noise level
Qnorm = Q %>%
  left_join(t0 %>% mutate(max0 = ifelse(max0 < noiseLevel, noiseLevel, max0))) %>%
  mutate(intensity = ifelse(digestTime > 0, intensity - max0, 0),
         intensity = ifelse(digestTime > 0 & intensity <= 0, NA, intensity))

remDataP = is.na(Qnorm$intensity) %>% which() %>% length()
paste0("remove ",remDataP, " out of ", nrow(Qnorm)," data points (", round(100*remDataP/nrow(Qnorm),2) ,"%) due to noise")

# data normalisation
# Q = Q %>%
#   group_by(pepSeq, biological_replicate, technical_replicate) %>%
#   mutate(intensity = intensity-intensity[digestTime == 0]) %>%
#   rowwise() %>%
#   mutate(intensity = ifelse(intensity < 0, 0, intensity)) %>%
#   ungroup()

# plotting
plotKinetics(Q, unlist(snakemake@output[["plot_raw"]]),
             meanTech = F, earlyOnly = F, sortByInt = T)
plotKinetics(Qnorm, unlist(snakemake@output[["plot_kinetics"]]),
             meanTech = T, earlyOnly = T, sortByInt = T)

# ----- generate output tables -----
# annotation
annotations = SKYLINE %>%
  group_by(pepSeq,substrateSeq) %>%
  summarise(noScans = n(),
            assignedScans = paste(ID,collapse = ";"),
            spectralAngles = paste(spectralAngle, collapse = ";"),
            qValues = paste(qValue,collapse = ";"),
            ionScores = paste(ionScore, collapse = ";"),
            charges = paste(charge, collapse = ";"),
            deltaRTs = paste(deltaRT/60, collapse = ";"),
            modifications = paste(modifiedSequence, collapse = ";"))

# kinetics
# mean over technical replicates
Qsum = Qnorm %>%
  group_by(pepSeq,biological_replicate,digestTime) %>%
  mutate(intensity_mean = mean(intensity, na.rm = T))

# summarise
kinetics = Qsum %>%
  filter(technical_replicate == 1) %>%
  group_by(substrateID,pepSeq,productType,spliceType,positions,biological_replicate) %>%
  summarise(digestTimes = paste(digestTime,collapse = ";"),
            intensities = paste(intensity_mean, collapse = ";"))

# join both
FINAL = full_join(annotations,kinetics) %>%
  select(substrateID,pepSeq,biological_replicate,digestTimes,intensities,
         substrateSeq,productType,spliceType,positions,noScans,assignedScans,
         spectralAngles,qValues,ionScores,deltaRTs,charges,modifications)


# statistics
FINAL$spliceType[FINAL$spliceType %in% c(NA,"")] = "PCP"
tbl = FINAL %>%
  uniquePeptides() %>%
  disentangleMultimappers.Type() %>%
  na.omit()

print(protein_name)

table(tbl$productType) %>%
  print()
table(tbl$spliceType) %>%
  print()


### OUTPUT ####
write.csv(FINAL, file = unlist(snakemake@output[["final_kinetics"]]),
          row.names = F)

Sys.time() %>% print()
