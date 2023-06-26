### aSPIRE ###
# description:  get peptide coverage map over time
# input:        TICs exported from Skyline, sample list
# output:       plotting of TICs
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(magick)
source("src/invitroSPI_utils.R")
source("src/aSPIre_plotting.R")
theme_set(theme_classic())

Mx = 10


print("-------------------------------------")
print("PLOTTING OF COVERAGE AND RESIDUE MAPS")
print("-------------------------------------")

### INPUT ###
Kinetics = read.csv(snakemake@input[["final_kinetics"]], stringsAsFactors = F)

### MAIN PART ###
# ----- preprocessing -----
# mean intensity over biological replicates
# assign weight to multi-mappers

DB = Kinetics %>%
  filter(!is.na(substrateID)) %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  dplyr::rename(digestTime = digestTimes,
                intensity = intensities) %>%
  dplyr::mutate(intensity = as.numeric(intensity),
                digestTime = as.numeric(digestTime)) %>%
  dplyr::group_by(substrateID, pepSeq, digestTime) %>%
  dplyr::mutate(intensity = mean(intensity)) %>%
  ILredundancy() %>%
  uniquePeptides() %>%
  resolve_multimapper() %>%
  tidyr::separate_rows(positions, sep = ";")


# ----- get coverage -----
getCoverage = function(protein_name) {
  
  print(protein_name)
  
  S = DB$substrateSeq[1]
  cntDB = DB
  
  no_pcp = cntDB$pepSeq[cntDB$productType == "PCP"] %>% unique() %>% length()
  no_psp = cntDB$pepSeq[cntDB$productType == "PSP"] %>% unique() %>% length()
  
  tp = cntDB$digestTime %>% unique() %>% as.numeric() %>% sort()
  tp = tp[which(tp != 0)]
  
  # --- get data
  # coverage for each time point
  D = lapply(tp, function(t){
    
    # print(t)
    
    cnt = cntDB[cntDB$digestTime == t, ]
    pos = str_split_fixed(cnt$positions, coll("_"), n=Inf)
    pos = apply(pos,2,as.numeric) %>%
      as.data.frame()
    
    d = data.frame(substrate = S,
                   residue = c(1:nchar(S)),
                   pcps = 0,
                   psps = 0,
                   n_pcps = 0,
                   n_psps = 0)
    
    for (i in d$residue) {
      
      d$n_pcps[i] = length(which(pos$V1 <= i & pos$V2 >= i & cnt$spliceType == "PCP"))
      d$n_psps[i] = length(which(pos$V1 <= i & pos$V2 >= i & cnt$spliceType != "PCP")) + length(which(pos$V3 <= i & pos$V4 >= i & cnt$spliceType != "PCP"))
      
      psp1 = cnt$intensity[which(pos$V1 <= i & pos$V2 >= i & cnt$spliceType != "PCP")] %>% log10() %>% sum(na.rm = T)
      psp2 = cnt$intensity[which(pos$V3 <= i & pos$V4 >= i & cnt$spliceType != "PCP")] %>% log10() %>% sum(na.rm = T)
      pcp = cnt$intensity[which(pos$V1 <= i & pos$V2 >= i & cnt$spliceType == "PCP")] %>% log10() %>% sum(na.rm = T)
      
      d$pcps[i] = pcp
      d$psps[i] = psp1+psp2
    }
    d[is.na(d)] = 0
    
    return(d)
  })
  
  names(D) = tp
  
  Dorig = plyr::ldply(D) %>%
    as.data.frame() %>%
    rename(digestTime = .id) %>%
    mutate(substrateID = protein_name)
  
  # normalise it
  Dorig$pcps = (Dorig$pcps-min(Dorig$pcps)) / (max(Dorig$pcps) - min(Dorig$pcps))
  Dorig$psps = (Dorig$psps-min(Dorig$psps)) / (max(Dorig$psps) - min(Dorig$psps))
  
  D = lapply(rev(tp), function(t){
    
    d = Dorig[Dorig$digestTime == t,]
    sp = smooth.spline(d$psps)
    nsp = smooth.spline(d$pcps)
    
    sp$y[sp$y < 0] = 0
    nsp$y[nsp$y < 0] = 0
    
    return(list(sp = sp, nsp = nsp))
  })
  names(D) = rev(tp)
  
  # --- plotting
  par(mfrow = c(2,2))
  
  lim = 1
  a = .3
  lw = 2
  xshift0 = -15
  yshift0 = -0.3
  L = nchar(S)
  
  # --- selected time point - spliced vs. non-spliced
  jsel = if(any(tp == 4)) which(tp == 4) else which.min(abs(4-tp))[1]
  
  selSP = D[[jsel]]$sp$y
  selSP = (selSP - min(selSP))/(max(selSP) - min(selSP))
  selNSP = D[[jsel]]$nsp$y
  selNSP = (selNSP - min(selNSP))/(max(selNSP) - min(selNSP))
  
  plot(x=D[[jsel]]$sp$x, y=-1*selSP,
       type = "l", lwd = lw,
       main = paste0(protein_name,"   - ",tp[jsel],"h"),
       sub = paste0("PCP: ",no_pcp,", PSP: ",no_psp),
       col = plottingCols["PSP"],
       ylim = c(-1.2*lim,1.2*lim),
       axes = F, xlab = "", ylab = "")
  polygon(c(D[[jsel]]$sp$x, rev(D[[jsel]]$sp$x)), c(-1*selSP, rep(0,length(D[[jsel]]$sp$y))),
          col = add.alpha(plottingCols["PSP"],a), lty = 0)
  # non-spliced
  lines(selNSP, col = plottingCols["PCP"], lwd = lw)
  polygon(c(D[[jsel]]$nsp$x, rev(D[[jsel]]$nsp$x)), c(selNSP, rep(0,length(D[[jsel]]$nsp$y))),
          col = add.alpha(plottingCols["PCP"],a), lty = 0)
  # legend and stuff
  abline(h = 0, lwd = 1)
  l = round(nchar(S),-1)
  axis(1, pos=-1.2*lim, at = seq(1, l+1, l/10), lwd = 0.5)
  axis(2, at = seq(-1,1,0.5), lwd = 0.5)
  
  
  plot.new()
  
  
  # --- spliced
  plot(x=D[[1]]$sp$x, y=D[[1]]$sp$y,
       type = "l", lwd = lw,
       main = "spliced", sub = "coverage over time",
       col = plottingCols["PSP"],
       ylim = c(yshift0*length(D),1.2*lim),
       xlim = c(xshift0*length(D), L+15),
       axes = F, xlab = "", ylab = "")
  polygon(c(D[[1]]$sp$x, rev(D[[1]]$sp$x)), c(D[[1]]$sp$y, rep(0,length(D[[1]]$sp$y))),
          col = add.alpha(plottingCols["PSP"],a), lty = 0)
  abline(h = 0, lwd = 1)
  
  xshift = xshift0
  yshift = yshift0
  for (j in 2:length(D)) {
    
    lines(x=D[[j]]$sp$x + xshift, y=D[[j]]$sp$y + yshift,
          col = plottingCols["PSP"], lwd = lw)
    polygon(c(D[[j]]$sp$x + xshift, rev(D[[j]]$sp$x + xshift)), c(D[[j]]$sp$y + yshift, rep(yshift,length(D[[j]]$sp$y))),
            col = add.alpha(plottingCols["PSP"],a), lty = 0)
    
    abline(h = 0+yshift, lwd = 1)
    
    xshift = xshift + xshift0
    yshift = yshift + yshift0
  }
  abline(b = yshift0/xshift0, a = -1*(yshift0/xshift0)*L)
  text(x = c(0, seq(xshift0, xshift-xshift0, xshift0))+L+10,
       y = c(0, seq(yshift0, yshift-yshift0, yshift0)),
       labels = rev(tp))
  
  
  # --- non-spliced
  # non-spliced
  plot(x=D[[1]]$nsp$x, y=D[[1]]$nsp$y,
       type = "l", lwd = lw,
       main = "non-spliced", sub = "coverage over time",
       col = plottingCols["PCP"],
       ylim = c(yshift0*length(D),1.2*lim),
       xlim = c(xshift0*length(D), L+15),
       axes = F, xlab = "", ylab = "")
  polygon(c(D[[1]]$nsp$x, rev(D[[1]]$nsp$x)), c(D[[1]]$nsp$y, rep(0,length(D[[1]]$nsp$y))),
          col = add.alpha(plottingCols["PCP"],a), lty = 0)
  abline(h = 0, lwd = 1)
  
  xshift = xshift0
  yshift = yshift0
  for (j in 2:length(D)) {
    
    lines(x=D[[j]]$nsp$x + xshift, y=D[[j]]$nsp$y + yshift,
          col = plottingCols["PCP"], lwd = lw)
    polygon(c(D[[j]]$nsp$x + xshift, rev(D[[j]]$nsp$x + xshift)), c(D[[j]]$nsp$y + yshift, rep(yshift,length(D[[j]]$nsp$y))),
            col = add.alpha(plottingCols["PCP"],a), lty = 0)
    
    abline(h = 0+yshift, lwd = 1)
    
    xshift = xshift + xshift0
    yshift = yshift + yshift0
  }
  
  abline(b = yshift0/xshift0, a = -1*(yshift0/xshift0)*L)
  text(x = c(0, seq(xshift0, xshift-xshift0, xshift0))+L+10,
       y = c(0, seq(yshift0, yshift-yshift0, yshift0)),
       labels = rev(tp))
  
  
  save(D, file = unlist(snakemake@output[["coveragevalues"]])) 
}


pdf(paste0("results/",protein_name,"/",protein_name,"_coverage.pdf"), height = 8, width = 12)
getCoverage(protein_name)
dev.off()
print("Coverage maps successfully plotted!....")


# ----- get residue maps -----
getResidueMap = function(protein_name) {
  
  print(protein_name)
  
  # --- get coverage
  load(unlist(snakemake@output[["coveragevalues"]]))
  tp = names(D) %>% as.numeric()
  jsel = if(any(tp == 4)) which(tp == 4) else which.min(abs(4-tp))[1]
  
  selSP = D[[jsel]]$sp$y
  selSP = (selSP - min(selSP))/(max(selSP) - min(selSP))
  selNSP = D[[jsel]]$nsp$y
  selNSP = (selNSP - min(selNSP))/(max(selNSP) - min(selNSP))
  
  selSP = Mx*selSP
  selNSP = Mx*selNSP
  
  # --- get peptides
  S = DB$substrateSeq[1]
  L = nchar(S)
  
  cntDB = DB %>%
    filter((digestTime == tp[jsel]))
  
  pos = str_split_fixed(cntDB$positions, coll("_"), n=Inf)
  pos = apply(pos,2,as.numeric) %>%
    as.data.frame()
  
  
  pdf(paste0("results/",protein_name,"/",protein_name,"_residuemap.pdf"), height = 4, width = 6)
  for (i in 1:nchar(S)) {
    
    k_pcp = which(pos$V1 <= i & pos$V2 >= i & cntDB$spliceType == "PCP")
    k_sr1 = which(pos$V1 <= i & pos$V2 >= i & cntDB$spliceType != "PCP")
    k_sr2 = which(pos$V3 <= i & pos$V4 >= i & cntDB$spliceType != "PCP")
    
    # get coordinates
    pcp = pos[k_pcp,c(1:2)]
    psp = pos[c(k_sr1, k_sr2),]
    v = c(rep("sr1", length(k_sr1)), rep("sr2", length(k_sr2)))
    
    # get intensities
    i_pcp = cntDB$intensity[k_pcp] %>% log10()
    i_psp = cntDB$intensity[c(k_sr1,k_sr2)] %>% log10()
    
    # --- plotting
    # normalise intensities
    i_pcp = selNSP[i]*(i_pcp/sum(i_pcp, na.rm = T))
    i_psp = selSP[i]*(i_psp/sum(i_psp, na.rm = T))
    
    # scaffold
    plot(x = 1:L, y = rep(NA,L), xlim = c(1,L),
         ylim = c(-Mx-2,Mx+2), axes = F, ylab = "peptides", xlab = "",
         main = protein_name, sub = paste0("residue ",i,", PCPs: ",nrow(pcp),", PSPs: ", nrow(psp)))
    abline(h = 0)
    
    # hotspot lines
    lines(x = 1:L, y = selNSP, lwd = 2, col = plottingCols["PCP"])
    lines(x = 1:L, y = -1*selSP, lwd = 2, col = plottingCols["PSP"])
    
    if (nrow(pcp) > 0 & nrow(psp) > 0) {
      
      # PCPs
      if (nrow(pcp) > 0) {
        ypos = 0
        for (j in 1:nrow(pcp)) {
          if (is.finite(i_pcp[j])) {
            ydiff = ypos+i_pcp[j]
            
            n = pcp[j,2]-pcp[j,1]+2
            xpos = (pcp[j,1]-0.5):(pcp[j,2]+0.5)
            
            polygon(x = c(xpos, rev(xpos)),
                    y = c(rep(ydiff,n), rep(ypos,n)),
                    col = plottingCols["PCP"], lty = 0)
            
            ypos = ydiff
          }
          
        } 
      }
      
      
      # PSPs
      if (nrow(psp) > 0){
        ypos = 0
        for (jj in 1:nrow(psp)) {
          if (is.finite(i_psp[jj])) {
            ydiff = ypos+i_psp[jj]
            
            n1 = psp[jj,2]-psp[jj,1]+2
            n2 = psp[jj,4]-psp[jj,3]+2
            
            xpos1 = (psp[jj,1]-0.5):(psp[jj,2]+0.5)
            polygon(x = c(xpos1, rev(xpos1)),
                    y = -1*c(rep(ydiff,n1), rep(ypos,n1)),
                    col = plottingCols["SR1"], lty = 0)
            
            xpos2 = (psp[jj,3]-0.5):(psp[jj,4]+0.5)
            polygon(x = c(xpos2, rev(xpos2)),
                    y = -1*c(rep(ydiff,n2), rep(ypos,n2)),
                    col = plottingCols["SR2"], lty = 0)
            
            ypos = ydiff
          }
          
        }
      }
    }
    
    # legend and stuff
    abline(v = i, lty = "dashed")
    axis(1)
    legend("topright",
           legend = c("PCP","SR1","SR2"), pch = rep(16,3),
           col = c(plottingCols["PCP"], plottingCols["SR1"], plottingCols["SR2"]),
           horiz = T, bty = "n", cex = .8)
    
  }
  dev.off()
  
}


getResidueMap(protein_name)
print("Residue maps successfully plotted!....")

print("Converting into GIF....")
# convert into gif
pdf_doc <- magick::image_read_pdf(paste0("results/",protein_name,"/",protein_name,"_residuemap.pdf"), density = 200)
gif_file_name_and_path <- unlist(snakemake@output[["residuemap"]])
pdf_doc %>% magick::image_animate(fps = 5) %>% magick::image_write(path = gif_file_name_and_path)
print("DONE!")


