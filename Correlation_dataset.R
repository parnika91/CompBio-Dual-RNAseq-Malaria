# This script performs permutation tests for gene coexpression analysis
# Running this once makes 10,000 permutations
# For each study, look for runs at different thresholds - all: all runs in the study, int: 50% transcriptome expression in host and in parasite, str: 70% transcriptome exp
# Put the datasets that can be analysed in Data/ 
# Make sure that there is disctinction made between datasets at different thresholds, such as with "<studyID>.ortho.data.int.RData"

library(WGCNA)
library(foreach)
library(doParallel)

# function to assign an RData object to a name
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

options(echo=TRUE)
args <- commandArgs()
studyID <- args[1] # for example, for a dataset called SRP1188996.ortho.data.int.RData, write SRP1188996.ortho.data.int in the command line
study <- loadRData(paste0("Data/", studyID, ".RData")) # if all the ortho.data.x files are stored in Data/ ; (x = int/str/all)
study <- t(study)

### correlation ####

# Original correlation coefficients, the test statistic
system.time(ori_cor <- cor(study, use = 'pairwise.complete.obs'))
n <- nrow(study)

reps <- 1000

PermAsso <- function()
{
  perm <- foreach(i=1:reps, .combine = "+") %do%
    {
      # Creating the null distribution here
      (abs(cor(study, study[sample(n, n),], use = 'pairwise.complete.obs') >= abs(ori_cor) +0))
    }
  return(perm)
}

PermAssoCompiled <- compiler::cmpfun(PermAsso)

outer_reps <- 10
registerDoParallel()

ptm <- proc.time()
outer <- foreach(j=1:outer_reps, .combine = "+") %dopar%
{ 
  print(j)
  PermAssoCompiled()
}
proc.time() - ptm

save(outer, file = paste0(studyID, "/outer_",studyID,"_",runif(1, min = 0, max = 10000),".RData", collapse = ''))
