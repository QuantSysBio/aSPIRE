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
print("subtract minimum intensity over all time points and set t=0 to 0")

Q = QUANTfiltered %>%
  tidyr::separate_rows(int_pasted, times_pasted, sep = ";") %>%
  rename(intensity = int_pasted,
         digestTime = times_pasted) %>%
  as.data.frame()

Q = Q %>%
  group_by(pepSeq, biological_replicate, technical_replicate) %>%
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(intensity = intensity-min(intensity))
Q$intensity[Q$digestTime == 0] = 0

plotKinetics(Q, unlist(snakemake@output[["plot_norm"]]),
             meanTech = F, earlyOnly = F, sortByInt = T)
plotKinetics(Q, unlist(snakemake@output[["plot_kinetics"]]),
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
Qsum = Q %>%
  group_by(pepSeq,biological_replicate,digestTime) %>%
  mutate(intensity_mean = if (all(intensity == 0) & all(digestTime != 0)) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T))

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
table(tbl$spliceType) %>%
  print()


### OUTPUT ####
write.csv(FINAL, file = unlist(snakemake@output[["final_kinetics"]]),
          row.names = F)

Sys.time() %>% print()
