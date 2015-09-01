# R script to create separate working directories for each OTU and put an OTU-specific data frame (containing all 
# M columns of metadata + 1 column of count data) into each working directory. 
# Also copy the R script for negative binomial model into each new working directory

# Clear workspace 
rm(list=ls())

# Load the data
df.leaf<-readRDS("df_leaf_common_rlog.rds")
df.root<-readRDS("df_root_common_rlog.rds")

# How many columns of metadata are there? default = 33
M<-33 # first M columns of the data frame are metadata
Y<-M+1  # Taxon counts will be stored in column number Y

### Make master tables of all OTUs to be tested in LEAVES and ROOTS, with integer keys ###
### The integer keys are to facilitate use of job_arrays on the cluster to run negative binomial models
leafKeys<-data.frame(OTU_ID=colnames(df.leaf)[Y:dim(df.leaf)[2]],leafKey=0)  # initialize with key = 0
rootKeys<-data.frame(OTU_ID=colnames(df.root)[Y:dim(df.root)[2]],rootKey=0)  # initialize with key = 0
leafKeys$leafKey<-rownames(leafKeys) # replace key with (integer) row names
rootKeys$rootKey<-rownames(rootKeys) # replace key with (integer) row names

# Make the 'leaf' and 'root' subdirectories
dir.create('./leaf')
dir.create('./root')
# Make directory to catch SLURM output files
dir.create('./slurm_outputs')

####### For each leaf taxon to be tested, make a subdirectory with an integer name (integers to enable easier processing as a Batch job array) and within the new directory, save a data frame that has all the necessary metadata and ONLY the relevant count data for the current taxon #######
for (i in 1:dim(leafKeys)[1]){
  dir.create(paste('./leaf/',i,sep='')) ## Create directory for current taxon
  thistaxon<-as.character(leafKeys[i,1]) ## Get name of current taxon
  thisdata<-df.leaf[,1:M] ## Extract metadata
  thisdata<-cbind(thisdata,df.leaf[[thistaxon]]) ## Add count data for current taxon only
  colnames(thisdata)[Y]<-thistaxon
  save(thisdata,file=paste('./leaf/',i,'/indata.RData',sep='')) ## Save reduced data frame in new directory (named as integer = job array index) for current taxon
  system(paste("cp countmodels_rlogLMM.R ./leaf/",i,"/countmodels_rlogLMM.R",sep="")) # Copy statistics script into new directory for current taxon
  rm(thistaxon,thisdata) # clean up
}

####### For each root taxon to be tested, make a directory with an integer name (integers to enable easier processing as a Batch job array) and within the new directory, save a data frame that has all the necessary metadata and ONLY the relevant count data for the current taxon #######
for (i in 1:dim(rootKeys)[1]){
  dir.create(paste('./root/',i,sep='')) ## Create directory for current taxon
  thistaxon<-as.character(rootKeys[i,1]) ## Get name of current taxon
  thisdata<-df.root[,1:M] ## Extract metadata
  thisdata<-cbind(thisdata,df.root[[thistaxon]]) ## Add count data for current taxon only
  colnames(thisdata)[Y]<-thistaxon # name the column for current taxon
  save(thisdata,file=paste('./root/',i,'/indata.RData',sep='')) ## Save reduced data frame in new directory for current taxon
  system(paste("cp countmodels_rlogLMM.R ./root/",i,"/countmodels_rlogLMM.R",sep="")) # Copy statistics script into new directory for current taxon
  rm(thistaxon,thisdata) # clean up
}

####### Write bash script to start job array for leaf taxa ######
sink('start_countmodels_rlogLMM_leaf.q')
cat("#!/bin/bash\n#\n")
cat(sprintf("#SBATCH --array=1-%d\n", dim(leafKeys)[1]))
cat("#SBATCH --job-name=leaf-%j \n")
cat("#SBATCH --output=./slurm_outputs/leaf-%j.slurmout \n")
cat("# Script to model individual taxa counts for Ecotypes microbiome assembly project\n")
cat("# See countmodels_rlogLMM.R for more information on statistical processes and outputs\n\n")
cat("cd ./leaf/$SLURM_ARRAY_TASK_ID \n") # change to integer working directory containing one taxon
cat("R CMD BATCH countmodels_rlogLMM.R \n") # execute statistics script in R
sink()

####### Write bash script to start job array for root taxa ######
sink('start_countmodels_rlogLMM_root.q')
cat("#!/bin/bash\n#\n")
cat(sprintf("#SBATCH --array=1-%d\n", dim(rootKeys)[1]))
cat("#SBATCH --job-name=root-%j \n")
cat("#SBATCH --output=./slurm_outputs/root-%j.slurmout \n")
cat("# Script to model individual taxa counts for Ecotypes microbiome assembly project\n")
cat("# See countmodels_rlogLMM.R for more information on statistical processes and outputs\n\n")
cat("cd ./root/$SLURM_ARRAY_TASK_ID \n") # change to integer working directory containing one taxon
cat("R CMD BATCH countmodels_rlogLMM.R \n") # execute statistics script in R
sink()

####### Submit bash scripts to SLURM queue #######
system("sbatch start_countmodels_rlogLMM_leaf.q")
system("sbatch start_countmodels_rlogLMM_root.q")
