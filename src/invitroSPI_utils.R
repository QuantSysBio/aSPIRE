### IN VITRO SPI - DOWNSTREAM ###
# description:  utils for DB analysis
# input:        -
# output:       -
# author:       HR

library(dplyr)
library(stringr)


########## plotting colours ########## 
plottingCols = c(
  PCP = "#EC9A56",
  PSP = "#7B80C7",
  cis = "darkslateblue",
  revCis = "lightskyblue",
  trans = "#BA69BE",
  allcis = "#9BBFE5",
  multimapper = "gray",
  
  SR1 = "#5C6095",
  SR2 = "#B0B3DD",
  
  PaesM = "gray",
  invitroSPI = "coral",
  qiSPI = "aquamarine",
  inSPIRE = "deeppink",
  
  proteins = "#405255",
  polypeptides = "#B9C997",
  
  SpechtDB = "olivedrab1",
  wholeDB = "olivedrab4",
  ProteasomeDB = "olivedrab4",
  randomDB = "gray80"
)



########## disentange rows containing multimappers ########## 
disentangleMultimappers.Type = function(DB, silent=T) {
  
  if(!silent) {
    print("DISENTANGLE MULTI-MAPPERS FOR PRODUCT TYPE")
  }
  
  
  # only considering PSP multi-mappers
  k = which(str_detect(DB$positions, coll(";")) & !str_detect(DB$productType, "PCP"))
  
  if (length(k) > 0) {
    
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    mm = list()
    
    # pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3, file = "tmp")
    for (r in 1:nrow(DB_mm)) {
      
      # setTxtProgressBar(pb, r)
      
      
      cnt_types = str_split(DB_mm$spliceType[r], coll(";"), simplify = T) %>%
        paste()
      
      cnt_pos = str_split(DB_mm$positions[r], coll(";"), simplify = T) %>%
        paste()
      cntDB = DB_mm[r, ]
      
      if (length(unlist(strsplit(cnt_pos[1],"_"))) == 2) {
        
        cntDB$productType = "PCP"
        cntDB$spliceType = "PCP"
        
      } else {
        
        # all multi-mappers correspond to the same product type
        if (cnt_types %>% unique() %>% length() == 1) {
          cntDB$spliceType = cnt_types %>% unique()
          
          # in case of cis and trans --> keep cis
        } else if ( any("trans" %in% cnt_types) ) {
          
          x = which(cnt_types == "trans")
          cnt_types = cnt_types[-x]
          cnt_pos = cnt_pos[-x]
          
          if (cnt_types %>% unique() %>% length() == 1) {
            cntDB$spliceType = cnt_types %>% unique()
            
          } else {
            cntDB$spliceType = "type_multi-mapper"
          }
          
          cntDB$positions = cnt_pos %>% paste(collapse = ";")
          
        } else {
          cntDB$spliceType = "type_multi-mapper"
          cntDB$positions = cnt_pos %>% paste(collapse = ";")
        }
        
      }
      
      mm[[r]] = cntDB
    }
    
    mm = plyr::ldply(mm)
    DB = rbind(DB_nomm, mm) %>%
      as.data.frame()
  }
  
  # beep("mario")
  return(DB)
}


disentangleMultimappers.AA = function(DB, retSinglePos = T) {
  
  print("DISENTANGLE MULTI-MAPPERS FOR AA AT sP1 AND sP1'")
  
  DB$AAmultimapper = "no"
  
  # only considering PSP multi-mappers
  k = which(str_detect(DB$positions, coll(";")) & !str_detect(DB$productType, "PCP"))
  
  if (length(k) > 0) {
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3)
    for (r in 1:nrow(DB_mm)) {
      
      setTxtProgressBar(pb, r)
      
      cnt_pos = strsplit(DB_mm$positions[r], ";") %>%
        unlist() %>%
        str_split_fixed(pattern = coll("_"), n = Inf)
      
      sP1 = str_sub(DB_mm$substrateSeq[r], start = cnt_pos[, 2], end = cnt_pos[, 2])
      sP1d = str_sub(DB_mm$substrateSeq[r], start = cnt_pos[, 3], end = cnt_pos[, 3])
      
      if ((length(unique(sP1)) > 1) | (length(unique(sP1d)) > 1)) {
        DB_mm$AAmultimapper[r] = "yes"
        
      } else if (retSinglePos) {
        DB_mm$positions[r] = paste(cnt_pos[1, c(1:4)], collapse = "_")
      }
      
    }
    
    DB = rbind(DB_nomm, DB_mm) %>%
      as.data.frame()
  }
  
  
  # beep("mario")
  return(DB)
}


disentangleMultimappers.SRlen = function(DB, retSinglePos = T) {
  
  print("DISENTANGLE MULTI-MAPPERS FOR SR length'")
  
  DB$SRmultimapper = "no"
  
  # only considering PSP multi-mappers
  k = which(str_detect(DB$positions, coll(";")) & !str_detect(DB$productType, "PCP"))
  
  if (length(k) > 0) {
    
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3)
    for (r in 1:nrow(DB_mm)) {
      
      setTxtProgressBar(pb, r)
      
      cnt_pos = strsplit(DB_mm$positions[r], ";") %>%
        unlist() %>%
        str_split_fixed(pattern = coll("_"), n = Inf)
      
      SR1len = as.numeric(cnt_pos[, 2]) - as.numeric(cnt_pos[, 1]) + 1
      SR2len = as.numeric(cnt_pos[, 4]) - as.numeric(cnt_pos[, 3]) + 1
      
      if ((length(unique(SR1len)) > 1) | (length(unique(SR2len)) > 1)) {
        DB_mm$SRmultimapper[r] = "yes"
        
      }  else if (retSinglePos) {
        DB_mm$positions[r] = paste(cnt_pos[1, c(1:4)], collapse = "_")
      }
      
    }
    
    DB = rbind(DB_nomm, DB_mm) %>%
      as.data.frame()
    
    
  }
  
  
  # beep("mario")
  return(DB)
}

disentangleMultimappers.IVSeqLen = function(DB, retSinglePos = T) {
  
  print("DISENTANGLE MULTI-MAPPERS FOR INTERVENING SEQUENCE length'")
  
  DB$IVseqmultimapper = "no"
  
  # only considering PSP multi-mappers
  k = which(str_detect(DB$positions, coll(";")) & !str_detect(DB$productType, "PCP"))
  
  if (length(k) > 0) {
    
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3)
    for (r in 1:nrow(DB_mm)) {
      
      setTxtProgressBar(pb, r)
      
      cnt_pos = strsplit(DB_mm$positions[r], ";") %>%
        unlist() %>%
        str_split_fixed(pattern = coll("_"), n = Inf)
      
      cnt_types = str_split(DB_mm$spliceType[r], coll(";"), simplify = T) %>%
        paste()
      
      ivLens = rep(NA, nrow(cnt_pos))
      for (p in 1:nrow(cnt_pos)) {
        
        # cis and trans
        if (cnt_pos[p,3] >= cnt_pos[p,1]) {
          ivLens[p] = (abs(as.numeric(cnt_pos[p, 3]) - as.numeric(cnt_pos[p, 2])) - 1) %>%
            as.numeric()
          
        # revCis
        } else {
          ivLens[p] = (abs(as.numeric(cnt_pos[p, 1]) - as.numeric(cnt_pos[p, 4])) - 1) %>%
            as.numeric()
        }
        
      }
      
      
      if (length(unique(ivLens)) > 1) {
        DB_mm$IVseqmultimapper[r] = "yes"
        
      }  else if (retSinglePos) {
        DB_mm$positions[r] = paste(cnt_pos[1, c(1:4)], collapse = "_")
      }
      
    }
    
    DB = rbind(DB_nomm, DB_mm) %>%
      as.data.frame()
    
  }
  
  # beep("mario")
  return(DB)
}


removeMultimappers.Type = function(DB) {
  print("REMOVE PEPTIDES THAT CAN SOLELY BE MULTI-MAPPERS FOR PRODUCT TYPE")
  
  k = which(DB$spliceType == "type_multi-mapper")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


removeMultimappers.AA = function(DB) {
  print("REMOVE PEPTIDES THAT ARE MULTI-MAPPERS IN TERMS OF AA AT sP1 OR sP1'")
  
  k = which(DB$AAmultimapper == "yes")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


removeMultimappers.SRlen = function(DB) {
  print("REMOVE PEPTIDES THAT ARE MULTI-MAPPERS IN TERMS OF SR LENGTH")
  
  k = which(DB$SRmultimapper == "yes")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

removeMultimappers.IVSeqLen = function(DB) {
  print("REMOVE PEPTIDES THAT ARE MULTI-MAPPERS IN TERMS OF INTERVENING SEQUENCE LENGTH")
  
  k = which(DB$IVseqmultimapper == "yes")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


########## filter for peptide length ########## 
filterPepLength = function(DB, cutoff=5) {
  print(paste0("REMOVE PEPTIDES THAT ARE SHORTER THAN ", cutoff, " AA"))
  
  k = which(nchar(DB$pepSeq) < cutoff)
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

########## synthesis errors ########## 

remSynthErrors = function(DB) {
  print("REMOVE SYNTHESIS ERRORS")
  
  k = which(str_detect(DB$productType, "synError"))
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


keepSynthErrors = function(DB) {
  print("FILTER FOR SYNTHESIS ERRORS")
  
  k = which(str_detect(DB$productType, "synError"))
  if (length(k) > 0) {
    DB = DB[k, ]
  }
  
  return(DB)
}

########## early time points only ########## 
filterEarlyTimepoints = function(DB) {
  print("FILTER FOR EARLY TIME POINTS")
  
  DB = DB[which(DB$digestTime %in% c(2, 4)), ]
  
  return(DB)
}


########## 20S standard proteasome only ########## 
filter20Sstandard = function(DB) {
  print("FILTER FOR 20S STANDARD PROTEASOME")
  
  if ("protIsotype" %in% names(DB)) {
    DB = DB[DB$protIsotype %in% c("20S standard", "20S K562"), ]
  }
  
  return(DB)
}


########## I/L redundancy in product sequences #########
ILredundancy = function(DB) {
  
  print("REPLACE ALL ISOLEUCINS BY LEUCINS")
  
  DB$pepSeq = str_replace_all(DB$pepSeq, "I", "L")
  DB$substrateSeq = str_replace_all(DB$substrateSeq, "I", "L")
  
  return(DB)
}

ILredundancy_pepsOnly = function(DB) {
  
  print("REPLACE ALL ISOLEUCINS BY LEUCINS")
  
  DB$pepSeq = str_replace_all(DB$pepSeq, "I", "L")
  
  return(DB)
}

########## remove peptides containing PTMS ########## 
remPTMs = function(DB) {
  print("REMOVE PTM CONTAINING PRODUCTS")
  
  k = which(! DB$PTM %in% c("", NA))
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

########## remove substrate sequence from PCPs ##########
removeSubstrateFromPCPs = function(DB) {
  
  k = which(str_detect(DB$productType, "PCP") & DB$substrateSeq == DB$pepSeq)
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

########## no products containing substrate N or C ########## 
remSubstrateNorCterm = function(DB) {
  
  print("REMOVE PRODCUCTS THAT CONTAIN THE SUBSTRATE'S N OR C TERM")
  
  rem = c()
  
  pb = txtProgressBar(min = 0, max = nrow(DB), style = 3)
  for (i in 1:nrow(DB)) {
    
    setTxtProgressBar(pb, i)
    
    pos = strsplit(DB$positions[i], ";") %>%
      unlist() %>%
      str_split_fixed(pattern = "_", Inf)
    
    substrateLen = nchar(DB$substrateSeq[i])
    
    cntNCterm = rep(F, nrow(pos))
    
    for (j in 1:nrow(pos)) {
      
      if (str_detect(DB$productType[i], "PCP")) {
        if ((pos[j, 1] == 1) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
      } else if (str_detect(DB$productType[i], "PSP")) {
        if ((pos[j, 1] == 1) | (pos[j, 3] == 1) | (pos[j, 4] == substrateLen) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
      }
      
    }
    
    if (all(cntNCterm)) {
      rem = c(rem, i)
    }
    
  }
  
  if (length(rem) > 0) {
    DB = DB[-rem, ]
  }
  
  return(DB)
}

########## no products containing SR of 1 aa length ########## 
remSR1aa = function(DB) {
  
  print("REMOVE PRODCUCTS THAT CONTAIN A SR OF 1 AA LENGTH")
  
  rem = c()
  
  # PSPs with 1 aa SR
  for (i in which(str_detect(DB$productType, "PSP"))) {
    
    cnt_pos = strsplit(DB$positions[i], ";") %>%
      unlist() %>%
      str_split_fixed(pattern = "_", Inf)
    
    cnt_sr1 = as.numeric(cnt_pos[, 2]) - as.numeric(cnt_pos[, 1]) + 1
    cnt_sr2 = as.numeric(cnt_pos[, 4]) - as.numeric(cnt_pos[, 3]) + 1
    
    if ( all(cnt_sr1 == 1) | all(cnt_sr2 == 1) ) {
      rem = c(rem, i)
    }
  }
  
  if (length(rem) > 0) {
    DB = DB[-rem, ]
  }
  
  return(DB)
}

########## label products containing substrate's N/C term or a SR of 1 aa length ##########

filterTypes = function(DB) {
  
  DB.spliceTypes = DB$spliceType
  
  pb = txtProgressBar(min = 0, max = nrow(DB), style = 3)
  
  for(i in 1:nrow(DB)) {
    
    setTxtProgressBar(pb, i)
    
    pos = strsplit(DB$positions[i], ";") %>%
      unlist() %>%
      str_split_fixed(pattern = "_", Inf)
    
    if (all(pos == "")) {
      pos = str_locate(DB$substrateSeq[i], DB$pepSeq[i])
    }
    
    substrateLen = nchar(DB$substrateSeq[i])
    
    # N/C term containing PCPs and PSPs
    cntNCterm = rep(F, nrow(pos))
    
    for (j in 1:nrow(pos)) {
      
      if (str_detect(DB$productType[i], "PCP")) {
        
        if ((pos[j, 1] == 1) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
        
      } else if (str_detect(DB$productType[i], "PSP")) {
        
        if ((pos[j, 1] == 1) | (pos[j, 3] == 1) | (pos[j, 4] == substrateLen) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
        
      }
      
    }
    
    if (all(cntNCterm)) {
      DB.spliceTypes[i] = paste0(DB.spliceTypes[i], "_NCterm")
    } else {
      DB.spliceTypes[i] = DB.spliceTypes[i]
    }
    
    # PSPs with 1 aa SR
    if (str_detect(DB$productType[i], "PSP")) { 
      
      cnt_sr1 = as.numeric(pos[, 2]) - as.numeric(pos[, 1]) + 1
      cnt_sr2 = as.numeric(pos[, 4]) - as.numeric(pos[, 3]) + 1
      
      if ( all(cnt_sr1 == 1) | all(cnt_sr2 == 1) ) {
        DB.spliceTypes[i] = paste0(DB.spliceTypes[i], "_1aaSR")
      } else {
        DB.spliceTypes[i] = DB.spliceTypes[i]
      }
      
    } else {
      DB.spliceTypes[i] = DB.spliceTypes[i]
    }
    
  }
  
  DB$spliceType = DB.spliceTypes %>% unlist()
  
  return(DB)
}


########## unique product sequences ########## 
uniquePeptides = function(DB) {
  print("UNIQUE PEPTIDES")
  
  DB = DB %>%
    distinct(substrateID, substrateSeq, pepSeq, productType,
             #spliceType,
             positions, .keep_all = T)
  
  # DB2 = DB %>%
  #   group_by(substrateID, substrateSeq, pepSeq, productType,
  #            spliceType, positions) %>% slice(1)
  
  # DB3 = DB[, c("substrateID", "substrateSeq", "pepSeq",
  #             "productType", "spliceType", "positions")] %>%
  #   unique()
  
  # DB4 = DB[!duplicated(DB[, c("substrateID", "substrateSeq", "pepSeq",
  #                             "productType", "spliceType", "positions")]),]
  
  return(DB)
}

uniquePeptides.perTP = function(DB) {
  print("UNIQUE PEPTIDES FOR EACH TIME POINT")
  
  DB = DB %>%
    distinct(substrateID, substrateSeq, digestTime, pepSeq, productType,
             spliceType, positions, .keep_all = T)
  
  
  return(DB)
}

########## cosmetics ########## 
DBcosmetics = function(DB) {
  
  print("COSMETICS AND STATS")
  
  if ("X" %in% names(DB)) {
    DB$X = NULL
  }
  
  # replace NA in splice type column of PCPs
  DB$spliceType[str_detect(DB$productType, "PCP")] = "PCP"
  
  # # statistics
  # print(paste0("number of identified products: ", nrow(DB)))
  # 
  # print(table(DB$spliceType))
  # print(table(DB$spliceType) / nrow(DB))
  # 
  # print("substrates")
  # print(paste0("number of substrates: ", DB$substrateSeq %>% unique() %>% length()))
  # DB$substrateID %>%
  #   unique() %>% 
  #   paste(collapse = ", ") %>%
  #   print()
  
  return(DB)
}

########## mapping functions ##########
getPositions <- function(seq,substrate){
  
  
  
  #########################
  # PCP
  #########################
  
  
  l = nchar(seq)
  
  k = which((grepl(pattern=seq,x=substrate)==TRUE))
  if(length(k)>0){
    
    pcp = numeric()
    
    for(j in 1:length(k)){
      a = substrate
      x = strsplit(a,split=seq)[[1]]
      nn = nchar(x)
      n1 = rep(NA,(length(nn)-1))
      n2 = rep(NA,(length(nn)-1))
      for(r in 1:(length(x)-1)){
        n1[r] = sum(nn[1:r])+(r-1)*nchar(seq)+1
        n2[r] = n1[r]+nchar(seq)-1
      }
      pcp = rbind(pcp,cbind(n1,n2))
    }
    return(pcp)
  }
  
  
  
  #########################
  # PSP
  #########################
  
  
  
  
  ll = nchar(seq)
  
  
  pept = unlist(seq)
  N = nchar(seq)
  
  # split peptides to P matrix
  P = strsplit(pept,split="")[[1]]
  
  # get permutations of length N
  x = c(1:N)
  y = c(1:N)
  z = as.vector(outer(x,y,paste,sep="_"))
  q = matrix(NA,length(z),2)
  for(i in 1:length(z)){
    q[i,] = as.numeric(strsplit(z[i],split="_")[[1]])
  }
  
  qs = apply(q,1,sum)
  k = which(qs==N)
  q = q[k,]
  
  # loop over all peptides
  res2 <- list()
  res1 <- list()
  
  psp <- list()
  
  psp <- list()
  res1 <- list()
  res2 <- list()
  
  # generate all strings for searches
  S = matrix(NA,dim(q)[1],2)
  for(i in 1:dim(q)[1]){
    S[i,1] = paste(P[1:q[i,1]],sep="",collapse="")
    S[i,2] = paste(P[(q[i,1]+1):N],sep="",collapse="")
  }
  
  # search each entry in prot for the two corresponding fragments and extract acc and positions
  
  for(i in 1:dim(S)[1]){
    
    psp[[i]] <- list()
    res1[[i]] = which((grepl(pattern=S[i,1],x=substrate)==TRUE))
    res2[[i]] = which((grepl(pattern=S[i,2],x=substrate)==TRUE))
    
    kk = which(res1[[i]]%in%res2[[i]])
    k = res1[[i]][kk]
    if(length(k)>0){
      
      for(j in 1:length(k)){
        
        a = substrate
        
        
        x = strsplit(a,split=S[i,1])[[1]]
        nn = nchar(x)
        n1 = rep(NA,(length(nn)-1))
        n2 = rep(NA,(length(nn)-1))
        for(r in 1:(length(x)-1)){
          n1[r] = sum(nn[1:r])+(r-1)*nchar(S[i,1])+1
          n2[r] = n1[r]+nchar(S[i,1])-1
        }
        #check if substrate Cterm==S[i,1]
        len = nchar(S[i,1])
        y = paste(strsplit(a,split="")[[1]][(nchar(a)-len+1):nchar(a)],collapse="")
        if(S[i,1]==y){
          n1 = c(n1,nchar(a)-len+1)
          n2 = c(n2,nchar(a))
        }
        tmp = unique(apply(cbind(n1,n2),1,paste,collapse="_"))
        tmp2 = matrix(as.numeric(unlist(strsplit(tmp,split="_"))),length(tmp),2,byrow=TRUE)
        n1 = tmp2[,1]
        n2 = tmp2[,2]
        
        x = strsplit(a,split=S[i,2])[[1]]
        nn = nchar(x)
        n3 = rep(NA,(length(nn)-1))
        n4 = rep(NA,(length(nn)-1))
        for(r in 1:(length(x)-1)){
          n3[r] = sum(nn[1:r])+(r-1)*nchar(S[i,2])+1
          n4[r] = n3[r]+nchar(S[i,2])-1
        }
        #check if substrate Cterm==S[i,2]
        len = nchar(S[i,2])
        y = paste(strsplit(a,split="")[[1]][(nchar(a)-len+1):nchar(a)],collapse="")
        if(S[i,2]==y){
          n3 = c(n3,nchar(a)-len+1)
          n4 = c(n4,nchar(a))
        }
        tmp = unique(apply(cbind(n3,n4),1,paste,collapse="_"))
        tmp2 = matrix(as.numeric(unlist(strsplit(tmp,split="_"))),length(tmp),2,byrow=TRUE)
        n3 = tmp2[,1]
        n4 = tmp2[,2]
        
        # get all internal combinations and keep only those with intv<=25
        
        z = as.vector(outer(n2,n3,paste,sep="_"))
        y = matrix(NA,length(z),2)
        for(zz in 1:length(z)){
          y[zz,] = as.numeric(strsplit(z[zz],split="_")[[1]])
        }
        intv = y[,2]-y[,1]-1
        x = which(intv<0)
        if(length(x)>0){ intv[x] = y[x,1]-y[x,2]+1-nchar(S[i,1])-nchar(S[i,2]) }
        x = which(intv<0)
        #    if(length(x)>0){ intv[x] = 1000 }
        
        select = which(intv<=5000)
        
        nnn = length(select)
        if(nnn>0){
          psp[[i]][[j]] = matrix(NA,nnn,5)
          
          for(j2 in 1:nnn){
            
            psp[[i]][[j]][j2,] = c(pept,y[select[j2],1]-nchar(S[i,1])+1,y[select[j2],1],y[select[j2],2],y[select[j2],2]+nchar(S[i,2])-1)
          }
        }
        
      }
      
    }
    
    
  }
  
  # unlist results and return as unique matrix with all possible explanations as rows
  x = unlist(psp)
  #  psp = matrix(x,length(x)/5,5,byrow=FALSE)
  
  res = numeric()
  for(i in 1:length(psp)){
    if(length(psp[[i]])>0){
      for(j in 1:length(psp[[i]])){
        res = rbind(res,psp[[i]][[j]])
      }
    }
  }
  
  
  # print(res)
  
  return(res)
  
  
}


# re-map peptides
mapping = function(DB) {
  
  d = DB %>%
    dplyr::select(substrateSeq, productType, spliceType, pepSeq) %>%
    dplyr::distinct()
  
  # pb = txtProgressBar(min = 0, max = dim(d)[1], style = 3, file="tmp")
  
  for(i in 1:dim(d)[1]){
    # setTxtProgressBar(pb, i)
    
    if(!(d$productType[i]=="CONT" | d$productType[i]=="CONT_synError")){
      s = gsub("I","L",as.vector(d$pepSeq[i]))
      substrate = gsub("I","L",as.vector(d$substrateSeq[i]))
      x = getPositions(s,substrate)
      
      
      #PCP
      if(dim(x)[2]==2){
        d$positions[i] = paste(apply(x,1,paste,collapse="_"),collapse=";")
      }
      
      
      #PSP
      if(dim(x)[2]>2){
        # print(i)
        if(dim(x)[2]==5 & dim(x)[1]>1){
          d$positions[i] = paste(apply(x[,-1],1,paste,collapse="_"),collapse=";")
        }
        if(dim(x)[2]==5 & dim(x)[1]==1){
          d$positions[i] = paste(x[,-1],collapse="_")
        }
        
        types = rep("cis",dim(x)[1])
        
        intv = as.numeric(x[,4])-as.numeric(x[,3])
        k = which(intv<=0)
        if(length(k)>0){
          types[k] = "revCis"
          
          k2 = which(as.numeric(x[k,2])<=as.numeric(x[k,5]) & as.numeric(x[k,3])>=as.numeric(x[k,4]))
          if(length(k2)>0){
            types[k[k2]] = "trans"
          }
          
        }
        
        d$spliceType[i] = paste(types,collapse=";")
      }
    }
    
  }
  
  DB$productType = NULL
  DB$spliceType = NULL
  DB$positions = NULL
  
  
  DB = left_join(DB, d) %>%
    as.data.frame()
  
  # re-assign the product type
  pos = strsplit(DB$positions, "_")
  wr = which(sapply(pos, length) == 2 & str_detect(DB$productType, "PSP"))
  if (length(wr) > 0) {
    DB$productType[wr] = str_replace(DB$productType[wr], "PSP", "PCP")
  }
  
  xr = which(sapply(pos, length) == 4 & str_detect(DB$productType, "PCP"))
  if (length(xr) > 0) {
    DB$productType[xr] = str_replace(DB$productType[xr], "PCP", "PSP")
  }
  
  return(DB)
}



########## load random DBs #########
loadRandom = function(fpath) {
  
  print("LOAD RANDOM DATABASES")
  
  # read cis, revCis, trans and PCP
  # concatenate and return list of data frames
  
  cis = read.csv(fpath[str_detect(fpath, "_cis")], stringsAsFactors = F)
  revCis = read.csv(fpath[str_detect(fpath, "_revCis")], stringsAsFactors = F)
  trans = read.csv(fpath[str_detect(fpath, "_trans")], stringsAsFactors = F)
  
  if (any(str_detect(fpath, "PCP_"))) {
    PCP = read.csv(fpath[str_detect(fpath, "_PCP")], stringsAsFactors = F)
  } else {
    PCP = NA
  }
  
  random = list(cis = cis,
                revCis = revCis,
                trans = trans,
                PCP = PCP)
  
  return(random)
}


# ----- resolve multi-mappers -----
# assign weight to position multi-mappers
resolve_multimapper = function(ProteasomeDB) {
  
  if (!"biological_replicate" %in% names(ProteasomeDB)) {
    ProteasomeDB$biological_replicate = 1
  }
  
  cnt = ProteasomeDB %>%
    tidyr::separate_rows(positions, sep=";") %>%
    dplyr::group_by(substrateID, pepSeq, digestTime, intensity, biological_replicate) %>%
    dplyr::summarise(n = n(),
                     positions = paste(positions, collapse = ";"),
                     intensity = intensity/n) %>%
    unique()
  
  DB = left_join(ProteasomeDB %>% dplyr::select(-intensity), cnt) %>%
    dplyr::select(-n)
  
  return(DB)
}

