### aSPIRE ###
# description:  plotting of TICs
# input:        TICs exported from Skyline, sample list
# output:       plotting of TICs
# author:       HPR

library(dplyr)
library(stringr)
library(data.table)

protein_name = snakemake@params[["protein_name"]]

print("---------------------")
print("CHROMATOGRAM PLOTTING")
print("---------------------")

### INPUT ###
sample_list = read.csv(snakemake@input[["sample_list"]], stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name == protein_name, ]

TIC = fread(snakemake@input[["tics"]]) %>%
  as.data.frame()


### MAIN PART ###
# get digestion times
ord = sample_list$digestTime %>% unique() %>% as.numeric() %>% sort()

TIC = left_join(TIC, sample_list %>% 
                  mutate(raw_file = gsub("\\.mgf","\\.raw",raw_file)) %>%
                  rename(FileName = raw_file))
TIC = TIC %>%
  rename(tp = digestTime,
         bioRep = biological_replicate) %>%
  group_by(protein_name,substrateID,tp,bioRep) %>%
  mutate(techRep = row_number()) %>%
  ungroup()


# ----- whole kinetics for each replicate -----
pdf(paste0("results/",protein_name,"/",protein_name,"_TICs.pdf"), height = 2*5, width = 2*6)
par(mfrow = c(2,2))

for (b in unique(TIC$bioRep)) {
  for (r in unique(TIC$techRep)) {
    
    cnt = TIC[TIC$bioRep == b & TIC$techRep == r, ]
    co = rainbow(n=nrow(cnt), alpha = .3)
    
    dt = lapply(1:nrow(cnt), function(i){
      data.frame(tp = strsplit(cnt$Times[i], split = ",") %>% unlist() %>% as.numeric(),
                 int = strsplit(cnt$Intensities[i], split = ",") %>% unlist() %>% as.numeric())
    })
    
    # sort
    names(dt) = cnt$tp
    dt = dt[na.omit(match(ord, names(dt)))]
    
    lim = sapply(dt,function(x){max(x$int)}) %>% max() %>% ceiling()
    
    # plot
    plot(x = dt[[1]]$tp, y = dt[[1]]$int,
         main = paste0(protein_name, ": rep ", b, " - inj ", r),
         ylim = c(0,lim),
         type = "l", lwd = 1, xlab = "RT [min]", ylab = "intensity", col = co[1])
    if (length(dt) > 1) {
      for (j in 2:length(dt)) {
        lines(dt[[j]]$tp, dt[[j]]$int, col = co[j])
      }
    }
    
    legend("topleft", legend = names(dt),
           col = co, lty = "solid", horiz = T, bty = "n", cex = .8)
    
    
  }
}


# ----- per time point for each replicate -----
for (t in ord) {
  
  cnt = TIC[TIC$tp == t, ]
  if (nrow(cnt) > 0) {
    co = hcl.colors(n=nrow(cnt)+1, palette = "viridis", alpha = .3)
    
    dt = lapply(1:nrow(cnt), function(i){
      data.frame(tp = strsplit(cnt$Times[i], split = ",") %>% unlist() %>% as.numeric(),
                 int = strsplit(cnt$Intensities[i], split = ",") %>% unlist() %>% as.numeric())
    })
    names(dt) = paste0(cnt$bioRep,"-",cnt$techRep)
    
    lim = sapply(dt,function(x){max(x$int)}) %>% max() %>% ceiling()
    
    # plot
    plot(x = dt[[1]]$tp, y = dt[[1]]$int,
         main = paste0(protein_name, ": ", t, " hours"),
         ylim = c(0,lim),
         type = "l", lwd = 1, xlab = "RT [min]", ylab = "intensity", col = co[1])
    
    if (length(dt) > 1) {
      for (j in 2:length(dt)) {
        lines(dt[[j]]$tp, dt[[j]]$int, col = co[j])
      }
    }
    
    legend("topleft", legend = names(dt),
           col = co[1:nrow(cnt)], lty = "solid", horiz = T, bty = "n", cex = .8)
    
  } else {
    plot.new()
  }
  
}

dev.off()


### OUTPUT ###
print("All chromatograms successfully plotted!....")
write("TICs plotted", file = unlist(snakemake@output[["chromatogram"]]))
