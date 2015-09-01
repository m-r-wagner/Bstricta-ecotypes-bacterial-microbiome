#!/usr/bin/Rscript

### Script to run relative log-transformed LMMs for a single taxon
### Working directory should include 1) this script, and 2) "indata.RData", which is a data frame (named 'thisdf') that contains M columns of sample metadata (including linear predictors), and a Yth column with count data for one taxon.
### One taxon per working directory.
### A job array batch script will submit this program to the queue FROM EACH working directory (working dirs have integer names)

# Clear workspace 
rm(list=ls())

library(lme4)
library(lmerTest)
library(lsmeans)

sessionInfo()

# This should be the same as in populate.R :
# How many columns of metadata are there? default = 19
M<-46 # first M columns of the data frame are metadata
Y<-M+1  # Taxon counts will be stored in column number Y

load('indata.RData',.GlobalEnv) # load data frame

print("Loading data...")
thisdf<-thisdata # store copy of thisdf
thistaxon<-colnames(thisdf)[Y] # store name of current taxon
save(thistaxon,file="taxon_name.RData") # for easy import

print("Defining functions...")
####### Function: r2.LMM() => estimate R^2 for a model as correlation between fitted and observed values #######
r2.LMM <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

print("Beginning linear mixed model of transformed count data...")
####### lmer: LMM for regularized log transformed count data #######
options(contrasts = c("contr.sum","contr.poly")) # contrasts sum to 0
lmer.model<-try(lmer((thisdf[[thistaxon]])~Genotype+Site+Genotype:Site+Harvested+Genotype*Age+(1|Genotype:Line)+(1|Site:Block)+(1|newPlate)+logObs,data=thisdf,REML=TRUE))

if (class(lmer.model)=="try-error"){
  sink("../model_errors.txt",append=TRUE)
  print(thistaxon)
  print("\n")
  sink()
}

LMM.stats<-as.data.frame(anova(lmer.model,type=3,ddf='Satterthwaite'))[,3:6]
LMM.stats<-cbind(data.frame("Taxon"=thistaxon, "Term"=rownames(LMM.stats)),LMM.stats) # add a column with the taxon name and term name
colnames(LMM.stats)<-c("Taxon","Term","NumDF","DenDF","F_or_Chisq","P_uncorrected")

# Likelihood ratio tests for random effects
lmer.model.noBlock<-try(lmer((thisdf[[thistaxon]])~Genotype+Site+Genotype:Site+Harvested+Genotype*Age+(1|Genotype:Line)+(1|newPlate)+logObs,data=thisdf,REML=FALSE))
lmer.model.noLine<-try(lmer((thisdf[[thistaxon]])~Genotype+Site+Genotype:Site+Harvested+Genotype*Age+(1|Site:Block)+(1|newPlate)+logObs,data=thisdf,REML=FALSE))
lmer.model.noPlate<-try(lmer((thisdf[[thistaxon]])~Genotype+Site+Genotype:Site+Harvested+Genotype*Age+(1|Genotype:Line)+(1|Site:Block)+logObs,data=thisdf,REML=FALSE))

LRT.Block<-as.data.frame(anova(lmer.model,lmer.model.noBlock))
LRT.Line<-as.data.frame(anova(lmer.model,lmer.model.noLine))
LRT.Plate<-as.data.frame(anova(lmer.model,lmer.model.noPlate))

Block.stats<-data.frame("Taxon"=thistaxon,"Term"="Block","NumDF"=LRT.Block$"Chi Df"[2],"DenDF"=Inf,"F_or_Chisq"=LRT.Block$Chisq[2], "P_uncorrected"=LRT.Block$P[2])
Line.stats<-data.frame("Taxon"=thistaxon,"Term"="Line","NumDF"=LRT.Line$"Chi Df"[2],"DenDF"=Inf,"F_or_Chisq"=LRT.Line$Chisq[2], "P_uncorrected"=LRT.Line$P[2])
Plate.stats<-data.frame("Taxon"=thistaxon,"Term"="Plate","NumDF"=LRT.Plate$"Chi Df"[2],"DenDF"=Inf,"F_or_Chisq"=LRT.Plate$Chisq[2], "P_uncorrected"=LRT.Plate$P[2])

LMM.stats<-rbind(LMM.stats,Block.stats,Line.stats,Plate.stats)
                 
write.table(LMM.stats,file=paste(thistaxon,"--LMM_stats.txt",sep=""),sep='\t',row.names=FALSE,col.names=TRUE) # save F test output
print ("LMM statistics have been written to file")

# print diagnostic plots
pdf(file=paste(thistaxon,"--LMM_residual_plot.pdf",sep=""))
plot(residuals(lmer.model)~fitted(lmer.model),col="dark grey") # plot residuals against fitted values
abline(h=0) # add horizontal line y=0
dev.off()
pdf(file=paste(thistaxon,"--LMM_QQnorm_plot.pdf",sep=""))
qqnorm(residuals(lmer.model))
qqline(residuals(lmer.model),distribution=qnorm)
dev.off()

print ("Diagnostic plots have been written to file")

####### Heritability estimates #######
# Fit model for each site
lmm.Jam<-try(lmer((subset(thisdf,Site=='Jam')[[thistaxon]])~Genotype+Harvested+Genotype*Age+(1|Genotype:Line)+(1|Block)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Jam'),REML=FALSE))
lmm.Mah<-try(lmer((subset(thisdf,Site=='Mah')[[thistaxon]])~Genotype+Harvested+Genotype*Age+(1|Genotype:Line)+(1|Block)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Mah'),REML=FALSE))
lmm.Sil<-try(lmer((subset(thisdf,Site=='Sil')[[thistaxon]])~Genotype+Harvested+Genotype*Age+(1|Genotype:Line)+(1|Block)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Sil'),REML=FALSE))
# Refit all models without genetic terms
lmm.noG.Jam<-try(lmer((subset(thisdf,Site=='Jam')[[thistaxon]])~Harvested+Age+(1|Block)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Jam'),REML=FALSE))
lmm.noG.Mah<-try(lmer((subset(thisdf,Site=='Mah')[[thistaxon]])~Harvested+Age+(1|Block)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Mah'),REML=FALSE))
lmm.noG.Sil<-try(lmer((subset(thisdf,Site=='Sil')[[thistaxon]])~Harvested+Age+(1|Block)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Sil'),REML=FALSE))

# Refit model across all sites, without genetic terms:
lmm.noG.all<-try(lmer((thisdf[[thistaxon]])~Site+Harvested+Age+(1|Site:Block)+(1|newPlate)+logObs,data=thisdf,REML=TRUE))


# calculate Heritability in each site as change in R^2 with vs without genetic information (Genotype, GxS, GxA, Line)
h2.Jam=r2.LMM(lmm.Jam)-r2.LMM(lmm.noG.Jam)
h2.Mah=r2.LMM(lmm.Mah)-r2.LMM(lmm.noG.Mah)
h2.Sil=r2.LMM(lmm.Sil)-r2.LMM(lmm.noG.Sil)
h2.all=r2.LMM(lmer.model)-r2.LMM(lmm.noG.all)

# initialize data frame
h2<-data.frame("Taxon"=character(),"Site"=character(),"H2"=numeric())

h2<-rbind(h2,
          data.frame("Taxon"=thistaxon,"Site"="Jam","H2"=h2.Jam),
          data.frame("Taxon"=thistaxon,"Site"="Mah","H2"=h2.Mah),
          data.frame("Taxon"=thistaxon,"Site"="Sil","H2"=h2.Sil),
          data.frame("Taxon"=thistaxon,"Site"="all","H2"=h2.all))

write.table(h2,file=paste(thistaxon,"--heritability.txt",sep=""),sep='\t',row.names=FALSE,col.names=TRUE) # save heritability estimates

print(" Heritability estimates have been written to file")
