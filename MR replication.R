
##Clean up R

rm(list=ls(all=TRUE))

####################################
##Install  packages
#install.packages("devtools")

library(devtools)

#devtools::install_github('MRCIEU/TwoSampleMR')
library(TwoSampleMR)

library(ggplot2)


getwd()

setwd("…/extra")


###################################

# BMI 

###################################


##First, read in file name with the non-standard headings

bmi_exp_dat<-read_exposure_data(
  filename="/...beta_snp_outcome_exp_bmi.csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file - 

bmi_exp_dat <- clump_data(bmi_exp_dat)

##


cost_out_dat<-read_outcome_data(
  filename="/...beta_snp_outcome_exp_bmi.csv",
  snps = bmi_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome",
  se_col = "se_outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data 

dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = cost_out_dat
)


######################## Main results ########################


res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
resivw<-mr(dat, method_list=c("mr_ivw"))


res.summary

bmi_het<-mr_heterogeneity(dat)

##################################################################################

########### Express two-sample MR results in terms of units of BMI rather than SD ########################
units.res.summary<- data.frame("outcome"=res.summary$outcome, "exposure"=res.summary$exposure, "method"=res.summary$method, "nsnp"=res.summary$nsnp, "b"=res.summary$b/4.6, "se"=res.summary$se/4.6, "pval"=res.summary$pval, "id.outcome"=res.summary$id.outcome, "id.exposure"=res.summary$id.exposure)      

units.res.summary

units10.res.summary<- data.frame("outcome"=res.summary$outcome, "exposure"=res.summary$exposure, "method"=res.summary$method, "nsnp"=res.summary$nsnp, "b"=res.summary$b*10/4.6, "se"=res.summary$se*10/4.6, "pval"=res.summary$pval, "id.outcome"=res.summary$id.outcome, "id.exposure"=res.summary$id.exposure)      

units10.res.summary

####################################



# SYSTOLIC BLOOD PRESSURE 
#USING International Consortium of Blood Pressure summary statistics
#(https://gwas.mrcieu.ac.uk/datasets/ieu-b-38/)

###################################



##First, read in my file name with the non-standard headings

sbp_exp_dat<-read_exposure_data(
  filename="/...beta_snp_outcome_exp_sbp.csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file - 

sbp_exp_dat <- clump_data(sbp_exp_dat)

##


cost_out_dat<-read_outcome_data(
  filename="/...beta_snp_outcome_exp_sbp.csv",
  snps = sbp_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome",
  se_col = "se_outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data 

dat2 <- harmonise_data(
  exposure_dat =sbp_exp_dat, 
  outcome_dat = cost_out_dat
)


######################## Main results for SBP ########################


res.summary2 <- mr(dat2, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
resivw<-mr(dat2, method_list=c("mr_ivw"))
sbp_het<-mr_heterogeneity(dat2)

res.summary2


########### Express two-sample MR results in terms of units of SD rather than mmhg ########################
sd.res.summary2<- data.frame("outcome"=res.summary2$outcome, "exposure"=res.summary2$exposure, "method"=res.summary2$method, "nsnp"=res.summary2$nsnp, "b"=res.summary2$b*20.7, "se"=res.summary2$se*20.7, "pval"=res.summary2$pval, "id.outcome"=res.summary2$id.outcome, "id.exposure"=res.summary2$id.exposure)      

sd.res.summary2

units10.res.summary2<- data.frame("outcome"=res.summary2$outcome, "exposure"=res.summary2$exposure, "method"=res.summary2$method, "nsnp"=res.summary2$nsnp, "b"=res.summary2$b*10, "se"=res.summary2$se*10, "pval"=res.summary2$pval, "id.outcome"=res.summary2$id.outcome, "id.exposure"=res.summary2$id.exposure)      

units10.res.summary2


###################################

# WAIST CIRCUMFERENCE 
#using the GIANT summary statistics (https://gwas.mrcieu.ac.uk/datasets/ieu-a-61/)
#(https://gwas.mrcieu.ac.uk/datasets/ieu-b-38/)

###################################



wc_exp_dat<-read_exposure_data(
  filename="/...beta_snp_outcome_exp_wc.csv",
  sep = ",",
  snp_col = "snp",
  beta_col = "betaexposure",
  se_col = "seexposure",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "pvalexposure",
  units_col = "unitsexposure",
  #gene_col = "Gene",
  samplesize_col = "samplesizeexposure"
)



##Clump this file - 

wc_exp_dat <- clump_data(wc_exp_dat)

##


cost_out_dat<-read_outcome_data(
  filename="/...beta_snp_outcome_exp_wc.csv",
  snps = wc_exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta_outcome",
  se_col = "se_outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p",
  units_col = "unit_output",
  #gene_col = "Gene",
  samplesize_col = "samplesize_out"
)


##Harmonise data 

dat3 <- harmonise_data(
  exposure_dat = wc_exp_dat, 
  outcome_dat = cost_out_dat
)



######################## Main results for waist circumference ########################


res.summary3 <- mr(dat3, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
wc_het<-mr_heterogeneity(dat3)
res.summary3




########### Express two-sample MR results in terms of units of WC rather than SD ########################

##comments - this might be in continuous units already - check 

units.res.summary3<- data.frame("outcome"=res.summary3$outcome, "exposure"=res.summary3$exposure, "method"=res.summary3$method, "nsnp"=res.summary3$nsnp, "b"=res.summary3$b/12.52, "se"=res.summary3$se/12.52, "pval"=res.summary3$pval, "id.outcome"=res.summary3$id.outcome, "id.exposure"=res.summary3$id.exposure)      

units.res.summary3

units10.res.summary3<- data.frame("outcome"=res.summary3$outcome, "exposure"=res.summary3$exposure, "method"=res.summary3$method, "nsnp"=res.summary3$nsnp, "b"=res.summary3$b*10/12.52, "se"=res.summary3$se*10/12.52, "pval"=res.summary3$pval, "id.outcome"=res.summary3$id.outcome, "id.exposure"=res.summary3$id.exposure)      

units10.res.summary3 

######################

#Summarize all results 

#######################################################################

#Waist circumference
res.summary3
units.res.summary3
units10.res.summary3 
wc_het


#BMI
res.summary
units.res.summary
units10.res.summary
bmi_het

#Systolic blood pressure
res.summary2
sd.res.summary2
units10.res.summary2
sbp_het




 
