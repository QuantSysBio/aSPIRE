### aSPIre ###
# description:  wait for Skyline results
# input:        -
# output:       -
# author:       HPR

library(dplyr)
library(R.utils)
numCPU = parallel::detectCores()

timeOut = 30
maxCrashes = 1

protein_name = snakemake@params[["protein_name"]]
sink(file = paste0("results/",protein_name,"/log.txt"), append = T, split = T)

print("-----------")
print("RUN SKYLINE")
print("-----------")


### INPUT ###
automatedSkyline = snakemake@params[["automatedSkyline"]]
raw_file_loc = snakemake@params[["raw_file_loc"]]

ssl = snakemake@input[["ssl"]]
fasta = snakemake@input[["fasta"]]
result_file = unlist(snakemake@output[["skyline"]])

blib = gsub(".ssl",".blib",basename(ssl))
Date = gsub(pattern="-",replacement="",x=as.character(Sys.Date()))


### MAIN PART ###
if (automatedSkyline) {
  
  # ----- AUTOMATED SKYLINE -----
  print("RUNNING AUTOMATED SKYLINE IN DOCKER")
  
  # system(paste0("cd ",raw_file_loc,"; find . -type f ! -name '*.raw' -exec rm -rf {} \;"))
  
  # ----- copy files
  # copy Skyline config and Skyline report
  cpw1 = file.copy(from = "data/Skyline/_Skyline-template.sky", to = raw_file_loc, overwrite = T)
  cpw2 = file.copy(from = "data/Skyline/_Skyline-report.skyr", to = raw_file_loc, overwrite = T)
  
  # ssl file and fasta file
  cpw3 = file.copy(from = ssl, to = raw_file_loc, overwrite = T)
  cpw4 = file.copy(from = fasta, to = raw_file_loc, overwrite = T)
  
  # library building command
  # lib_command = paste0("docker run -i --rm -v ",raw_file_loc,":/data ",
  #                      "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses ",
  #                      "wine blibbuild ", basename(ssl)," ", blib)
  
  # Skyline command
  skyline_command = paste0("docker run -i --rm -v ",raw_file_loc,":/data ",
                           "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses ",
                           "wine SkylineCmd --timestamp --dir=/data --in=_Skyline-template.sky ",
                           "--save ",
                           "--out=",protein_name,"_",Date,".sky ",
                           
                           # build library and import results
                           # "--add-library-path=",blib," ",
                           # "--add-library-name=",gsub(".blib","",blib)," ",
                           # "--import-all-files=./ ",
                           # "--import-filename-pattern=raw ",
                           "--import-search-file=",basename(ssl)," ",
                           "--import-search-add-mods ",
                           "--import-search-include-ambiguous ",
                           
                           # load fasta
                           "--import-fasta=",basename(fasta)," ",
                           "--keep-empty-proteins ",
                           
                           # load transition list (results) - extracting chromatograms
                           # "--exp-file=",basename(ssl)," ",
                           # "--exp-strategy=single ",
                           # "--exp-method-type=standard ",
                           # "--exp-max-trans=1000000 ",
                           
                           # document settings
                           "--import-threads=",numCPU," ",
                           "--refine-auto-select-peptides ",
                           "--refine-auto-select-transitions ",
                           "--refine-auto-select-precursors ",
                           "--refine-min-peptides=0 ",
                           "--refine-min-transitions=0 ",
                           "--refine-min-peak-found-ratio=0 ",
                           "--refine-max-peak-found-ratio=1 ",
                           "--refine-minimum-detections=0 ",
                           "--refine-min-dotp=0 ",
                           "--refine-min-idotp=0 ",
                           "--tran-precursor-ion-charges=1,2,3,4,5,6 ",
                           "--tran-product-ion-types=y,p ",
                           "--full-scan-rt-filter-tolerance=0.5 ",
                           
                           # export TICs and report
                           "--chromatogram-file=",protein_name,"_TICs.tsv ",
                           "--chromatogram-tics ",
                           "--report-add=_Skyline-report.skyr --report-name=MS1_HPR --report-invariant ",
                           "--report-file=MS1_HPR.csv > ", dirname(ssl),"/skyline_log.txt")
  
  # ----- run Skyline in Docker
  print("running Skyline (this can take up to an hour) ....")
  if (cpw1 & cpw2 & cpw3 & cpw4) {
    # print(lib_command)
    # system(lib_command)
    
    print(skyline_command)
    
    no_crashes = 0
    while (no_crashes <= maxCrashes) {
      
      # run Skyline
      withTimeout(
        expr = system(skyline_command),
        timeout = timeOut*60
      )
      
      # copy results
      cpw5 = file.copy(from = paste0(raw_file_loc,"/",basename(result_file)), to = dirname(result_file), overwrite = T)
      cpw6 = file.copy(from = paste0(raw_file_loc,"/",protein_name,"_TICs.tsv"), to = dirname(result_file), overwrite = T)
      cpw7 = file.copy(from = paste0(raw_file_loc,"/",protein_name,"_",Date,".sky"), to = dirname(result_file), overwrite = T)
      
      # decide whether to try again
      if (!(cpw5&cpw6&cpw7)) {
        no_crashes = no_crashes+1
      } else {
        break
      }
      
    }
    
    print("Finished Skyline!")
  } else {
    print("Cannot run Skyline")
  }
  
  if (cpw5 & cpw6 & cpw7) {
    print("Success! All data are in place to continue")
  } else {
    print("WARNING! Results files could not be copied back, check your server connection or folder permission!")
    break
  }
  
  
} else {
  
  # ----- MANUAL SKYLINE -----
  print("RUNNING SKYLINE MANUALLY")
  
  # ----- user info
  print("------------------------------------------------------------------")
  print("please copy the following files to your .raw file location:")
  print(ssl)
  print(fasta)
  paste0("run Skyline and copy the following report to: results/",protein_name,": ")
  print(basename(result_file))
  print("------------------------------------------------------------------")
  
  # ----- wait -----
  while(!file.exists(result_file)) {
    
    paste0(Sys.time()," - waiting for Skyline output..") %>% print()
    Sys.sleep(30)
    
  }
  
}
