### aSPIre ###
# description:  plotting function for kinetics
# input:        -
# output:       -
# author:       HPR

library(dplyr, warn.conflicts = F)
library(stringr)
library(RColorBrewer)
library(wesanderson)
source("src/invitroSPI_utils.R")
options(dplyr.summarise.inform = FALSE)

requ = c("pepSeq","spliceType","positions","biological_replicate")


plotKinetics = function(Qtable, outfile, meanTech=F, earlyOnly=T, sortByInt=T) {
  
  # ----- preprocessing -----
  if (!all(requ %in% names(Qtable))) {
    print("!!! plotting requires the following information:")
    print(requ)
  }
  
  Qtable = disentangleMultimappers.Type(Qtable)
  Qtable$spliceType[is.na(Qtable$spliceType)] = "PCP"
  
  if ("int_pasted" %in% names(Qtable)) {
    Qtable = Qtable %>%
      tidyr::separate_rows(int_pasted, times_pasted, sep = ";") %>%
      rename(intensity = int_pasted,
             digestTime = times_pasted) %>%
      as.data.frame()
  }
  
  # ----- sort by intensity -----
  if (sortByInt) {
    x = Qtable %>%
      filter(digestTime == max(digestTime)) %>%
      arrange(desc(intensity))
    pepU = unique(x$pepSeq)
    
  } else {
    pepU = Qtable$pepSeq %>% unique()
  }
  
  pdf(outfile, height = 16, width = 16)
  par(mfrow = c(4,4))
  
  # ----- iterate peptides -----
  for (i in 1:length(pepU)) {
   
    cnt = Qtable[Qtable$pepSeq == pepU[i],] 
    cnt$intensity = as.numeric(cnt$intensity)
    cnt$digestTime = as.numeric(cnt$digestTime)
    
    # ----- mean over technical replicates -----
    if(meanTech) {
      tmp = cnt %>%
        group_by(pepSeq,spliceType,positions,biological_replicate,digestTime) %>%
        summarise(mean_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T),
                  sd_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else sd(intensity[intensity!=0 | digestTime == 0], na.rm=T)) %>%
        arrange(digestTime) %>%
        suppressMessages()
    } else {
      tmp = cnt %>%
        select(pepSeq,spliceType,positions,biological_replicate,digestTime,intensity) %>%
        mutate(mean_int = intensity,
               sd_int = 0) %>%
        arrange(digestTime) %>%
        suppressMessages()
    }
    k = which(!is.na(tmp$mean_int))
    
    if (earlyOnly) {
      tmp = tmp[tmp$digestTime < 20,]
    }
    
    # ----- add colours -----
    nCol = tmp %>% group_by(biological_replicate,digestTime) %>% summarise(n = n())
    nCol = nCol$n[1]*length(unique(tmp$biological_replicate))
    if (nCol <= 4) {
      Cols = wes_palette("GrandBudapest1",n=nCol,type = "discrete")
    } else {
      Cols = c(wes_palette("GrandBudapest1",n=4,type = "discrete"),
               wes_palette("GrandBudapest2",n=4,type = "discrete"),
               wes_palette("Moonrise2",n=4,type = "discrete"),
               wes_palette("IsleofDogs1",n=4,type = "discrete"),
               wes_palette("Royal1",n=4,type = "discrete"))
      Cols = Cols[1:nCol]
    }
    tmp = tmp %>%
      group_by(biological_replicate,digestTime) %>%
      mutate(technical_replicate = row_number(),
             combo = paste0(biological_replicate,"-",technical_replicate))
    tmp = data.frame(combo = unique(tmp$combo), col = as.character(Cols)) %>%
      right_join(tmp) %>%
      suppressMessages()
    
    # ----- dot plot -----
    plot(x=tmp$digestTime, y=tmp$mean_int,
         pch= 16, cex = 1.5,
         xlab = "time [hrs]", ylab = "intensity",
         col = tmp$col,
         main = tmp$positions[1],
         sub = paste0(tmp$spliceType[1], " - ", pepU[i]))
    
    # ----- add lines -----
    for (gr in unique(tmp$combo)) {
      lines(x=tmp$digestTime[tmp$combo == gr],
            y=tmp$mean_int[tmp$combo == gr],
            col = tmp$col[tmp$combo==gr][1], lwd=1.5)
      
      # add standard dev
      if(meanTech) {
        arrows(tmp$digestTime[tmp$combo == gr], tmp$mean_int[tmp$combo == gr]-tmp$sd_int[tmp$combo == gr],
               tmp$digestTime[tmp$combo == gr], tmp$mean_int[tmp$combo == gr]+tmp$sd_int[tmp$combo == gr],
               length=0.03, angle=90, code=3,
               lty = "solid", lwd = 1, col = tmp$col[tmp$combo==gr][1]) %>%
          suppressWarnings()
      }
      
    }
    
    # ----- add legend -----
    txt = unique(tmp$combo)
    if(meanTech) {txt = str_remove_all(txt,"-[:digit:]")}
    legend("topleft",
           legend = txt, lty = rep("solid",length(txt)), col = Cols,bty="n")
    
  }

  dev.off()
  
}

