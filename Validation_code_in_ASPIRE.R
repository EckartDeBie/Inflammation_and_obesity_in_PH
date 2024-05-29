#Sheffield data validation
library(googledrive)
library(googlesheets4)
library(readxl)
library(tidyverse)
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
googledrive::drive_auth()
#drive_find()
dir = drive_find(pattern = 'aspire_shared', type='folder')
query = paste('"', dir$id, '"',  ' in parents', sep='')
a <- drive_find(q=query)
#file needs to be google sheet! 
meta <- gs4_get("drive_location1")
downloaded_file <- read_sheet(meta)

#this works!
df <- downloaded_file
df2 <- df %>% separate_longer_delim(c(CRP), delim= "{n}")
df2[c('measurement', 'date', 'level', 'unit')] <- str_split_fixed(df2$CRP, ' ', 4)

#merge with clinical data
meta2 <- gs4_get("drive_location2")
treatment_response <- read_sheet(meta2)

#merge on nearest data
df2$date <- as.Date(df2$date, format = "%d/%m/%Y")
treatment_response$DateofFinalDiagnosis <- as.Date(treatment_response$DateofFinalDiagnosis, format='%Y-$m-%d')
df2$ID <- df2$HOSPITAL_ID
treatment_response$ID <- treatment_response$STHNumber

file_merge <- merge(df2, treatment_response, by.x = 'ID', by.y='STHNumber')
df_done <- file_merge

#rename some variables (only want numeric vars)
df_done$level <- ifelse(df_done$level %in% c('SB', 'RNA', 'BIS'), NA, df_done$level)
df_done$level <- ifelse(df_done$level == '>700', 700, df_done$level)
df_done$level <- ifelse(df_done$level == '<2', 2, df_done$level)
df_done$level <- ifelse(df_done$level == '<0.3', 0.3, df_done$level)
#filter out missing CRPs
df_done <- df_done %>% filter(!is.na(date))
#str trim
df_done$level <- str_trim(df_done$level)

#make numeric
df_done$level<- as.numeric(df_done$level)

#now remove NAs
df_done <- df_done %>% filter(!is.na(level))

#now get CRP closest to diagnostic date
df_done$time_diff <- abs(df_done$date - df_done$DateofFinalDiagnosis)
df_done$num_time_diff <- as.numeric(df_done$time_diff)

df_done_filt <- df_done %>% group_by(ID) %>% slice(which.min(num_time_diff))
crp_df <- df_done_filt

meta4 <- gs4_get("drive_location3")
mri_correct <- read_sheet(meta4)

#=====================================================================================
#data is all set up and ready to use
  
#load relevant packages
library(tidyverse)
library(readxl)
library(ggplot2)
library(survival)
library(survminer)


#===========================================================================
#start with CRP
#we need: id, type of PH, RA area, RAP, mPAP, PVR, PAWP, CO, age at diagnosis, sex, TAPSE, event (for survival), survival time subtype of PH, 6MWD, CRP, BMI, if possible: comorbidities
df_done2 <- df_done %>% select(ID, date, level, Gender, DateofFinalDiagnosis, AgeAtDiag, FinalPrimaryPHDiagnosis, LifeStatus, YearsDiagDeath, BMI, mRAP, mPAP, PAWP, PVR, CardiacOutput, FVCp, eGFR, BaseISWD)
df_done3 <- mri_correct %>% filter(mri_order ==1) %>% select(sth_id, who_functional_class, walking_distance)


df_crp <- merge(df_done2, df_done3, by.x='ID', by.y='sth_id', all=T)
#NOTE: slight differences between BaseISWD and walking distance exist --> be mindful of these differences! 

#==============================================================================
#STEP 2
#=============================================================================
#perform linear modelling of variables
#plot variables first to check if they make sense

#get a model with just diagnostic CRP data - as this is what we'll need for the model
diag_crp <- merge(crp_df[,c(1,5)], df_crp, by =c('ID','date'))
diag_crp <- unique(diag_crp)

#some IDs seem merged....
n_occur <- data.frame(table(diag_crp$ID))
ids_duplicated <- n_occur[n_occur$Freq > 1,]
#these people have 2 CRPs per day --> so take avarage CRP value for the day
df_double <- df_crp %>% filter(ID %in% ids_duplicated$Var1)
df_double <- merge(crp_df[,c(1,5)], df_double, by =c('ID','date'))
df_double2 <- df_double %>% group_by(ID) %>% summarise(mean_level = mean(level))
df_double2 <- merge(df_double, df_double2)
df_double2$level <- df_double2$mean_level
df_double2 <- df_double2[,-c(21)]
df_double2 <- unique(df_double2)

#now merge this wih df_crp
diag_crp <- merge(crp_df[,c(1,5)], df_crp, by =c('ID','date'))
diag_crp <- diag_crp %>% filter(!ID %in% ids_duplicated$Var1)
diag_crp <- rbind(diag_crp, df_double2)
diag_crp <- unique(diag_crp)
#now the number is correct! 


summary(as.factor(diag_crp$FinalPrimaryPHDiagnosis))

#filter for Lm on just adults and no PVOD due to low (n=2) numbers
clean_df <- diag_crp %>% filter(AgeAtDiag >=18)
clean_df <- clean_df %>% filter(FinalPrimaryPHDiagnosis != 'C - Pulmonary Veno Occlusive Disease and/or Pulmonary Capillary Haemangiomatosis')

#now investigate the distributions
hist(clean_df$level)
hist(clean_df$PVR)
hist(clean_df$BMI)
hist(clean_df$mRAP)
hist(clean_df$mPAP)
hist(clean_df$CardiacOutput)
hist(clean_df$BaseISWD)
hist(clean_df$AgeAtDiag)
hist(clean_df$FVCp)
hist(clean_df$eGFR)
hist(clean_df$BaseISWD)
hist(clean_df$PAWP)

#now generate a function to do linear modelling
#picking model which was optimal in IPAH for this

#now fit some LMs
library(ggResidpanel)

clean_df$Gender <- as.factor(clean_df$Gender)
clean_df$FinalPrimaryPHDiagnosis <- as.factor(clean_df$FinalPrimaryPHDiagnosis)

#also eGFR of 0 is not possible
clean_df$eGFR[clean_df$eGFR <= 0] <- NA

df <- clean_df

library(broom)

#write a formula
#NOTE HERE: WE NEED TO DECIDE ON LOG TRANSFORMING
do_lm_log = function(myvar) {
  
  #1) now model myvar
  #only take complete data to model the same df
  df2 <- df %>% select(ID, myvar, AgeAtDiag, Gender, BMI, level)
  df2 <- na.omit(df2)
  
  fit_3 <- lm(log(df2[[myvar]] +1)  ~ df2$AgeAtDiag + df2$Gender + log(df2$level) + log(df2$BMI))
  plot(resid_panel(fit_3,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  print(summary(fit_3))
  
  output1 <- tidy(summary(fit_3))
  output1$variable_tested <- myvar
  return(output1)
  
}

do_lm_no_log = function(myvar) {
  
  #1) now model myvar
  #only take complete data to model the same df
  df2 <- df %>% select(ID, myvar, AgeAtDiag, Gender, BMI, level)
  df2 <- na.omit(df2)
  
  fit_3 <- lm(df2[[myvar]]  ~ df2$AgeAtDiag + df2$Gender + log(df2$level) + log(df2$BMI))
  plot(resid_panel(fit_3,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  print(summary(fit_3))
  
  output1 <- tidy(summary(fit_3))
  output1$variable_tested <- myvar
  return(output1)
  
}

#to log-transform: PVR, eGFR, ISWD, PAWP
#add in CO, CO, PAP, NT-pro-BNP --> NT-pro-BNP not in ASPIRE
pvri <- df_done_filt %>% select(ID, date, PAWP, mPAP, CardiacIndex, CardiacOutput)
pvri$pvri <- (pvri$mPAP - pvri$PAWP)/pvri$CardiacIndex
pvri$tpr <- pvri$mPAP/pvri$CardiacOutput

#load in new data
meta5 <- gs4_get("drive_location4")
pyhon_sheet <- read_sheet(meta5)
bnp <- pyhon_sheet %>% select(sth_id, date_diagnosis, bnp_date, bnp)
bnp <- na.omit(bnp)
#select BNP close to diagnosis
bnp$time_diff <- abs(bnp$bnp_date - bnp$date_diagnosis)
#time diff is recorded in seconds
bnp$time_diff <- as.numeric(bnp$time_diff)
#take BNPs within a week of diagnosis
bnp <- bnp %>% filter(time_diff <= 604800)
#sometimes multiple measurements in a day --> take mean
bnp <- bnp %>% group_by(sth_id) %>% summarise(bnp_level = mean(bnp))

ci <- pvri %>% select(ID, CardiacIndex)
#merge data together
ci <- merge(ci, bnp, by.x='ID', by.y='sth_id', all=T)

df <- merge(df, ci, by='ID')
#not to log-transform: FVC, RAP
#plot distributions of other vars
hist(log(df$bnp_level))
#looks better with log
hist(log(df$CardiacIndex))
#looks better with log
hist(log(df$CardiacOutput))
#looks better with log
hist(log(df$mPAP))

var_nolog <- c('FVCp', 'mRAP', 'eGFR')
var_log <- c('PVR', 'BaseISWD', 'PAWP', 'mPAP','CardiacOutput', 'CardiacIndex', 'bnp_level')
#note, do add +1 for log-transformation

ao <- lapply(var_nolog, do_lm_no_log)
ao2 <- do.call(rbind, ao)

ab <- lapply(var_log, do_lm_log)
ab2 <- do.call(rbind, ab)

lm_1_output <- rbind(ao2, ab2)


#do ISWD model per group
df1 <- df %>% filter(FinalPrimaryPHDiagnosis == 'B - Pulmonary Arterial Hypertension')
df2 <- df %>% filter(FinalPrimaryPHDiagnosis == 'D - Pulmonary Hypertension due to Left Heart Disease')
df3 <- df %>% filter(FinalPrimaryPHDiagnosis == 'E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia')
df4 <- df %>% filter(FinalPrimaryPHDiagnosis == 'F - Chronic Thromboembolic Pulmonary Hypertension')
df5 <- df %>% filter(FinalPrimaryPHDiagnosis == 'G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms')

#function model for df
do_lm_per_group_iswd = function(df) { 
fit <- lm(log(df$BaseISWD + 10)  ~ df$AgeAtDiag + df$Gender + log(df$level) + log(df$BMI))
plot(resid_panel(fit,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
print(summary(fit)) 
}
do_lm_per_group_iswd(df1)
do_lm_per_group_iswd(df2)
do_lm_per_group_iswd(df3)
do_lm_per_group_iswd(df4)
do_lm_per_group_iswd(df5)
#================================================================
#STEP 3
#================================================================
#model BMI and CRP separately

#also model BMI (different set of analysis so no function needed)
bmi <- df %>% select(ID, BMI, level, AgeAtDiag, Gender)
bmi <- na.omit(bmi)

fit_bmi2 <- lm(log(bmi$BMI) ~ bmi$AgeAtDiag + as.factor(bmi$Gender) + log(bmi$level))
a1 <- tidy(summary(fit_bmi2))
a1$variable_tested <- 'BMI'
resid_panel(fit_bmi2,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE)

#also model crp (likewise as for BMI)
fit_crp2 <- lm(log(bmi$level) ~ bmi$AgeAtDiag + as.factor(bmi$Gender) + log(bmi$BMI))
a2 <- tidy(summary(fit_crp2))
a2$variable_tested <- 'CRP'
resid_panel(fit_crp2,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE)

abc123 <- rbind(a1,a2)

lm_output <- rbind(abc123, lm_1_output)
saveRDS(lm_output, 'summary_linear_models_ASPIRE.rds')


lm_res_adj <- lm_output
sd_crp <- sd(log(df$level +1))
sd_bmi <- df %>% filter(!is.na(BMI))
sd_bmi <- sd(log(sd_bmi$BMI +1))
sd_pawp <- df %>% filter(!is.na(PAWP))
sd_pawp <- sd(log(sd_pawp$PAWP + 1))
sd_pvr <- df %>% filter(!is.na(PVR))
sd_pvr <- sd(log(sd_pvr$PVR + 1))
sd_iswd <- df %>% filter(!is.na(BaseISWD))
sd_iswd <- sd(log(sd_iswd$BaseISWD + 1))
sd_co <- df %>% filter(!is.na(CardiacOutput))
sd_co <- sd(log(sd_co$CardiacOutput + 1))
sd_ci <- df %>% filter(!is.na(CardiacIndex))
sd_ci <- sd(log(sd_ci$CardiacIndex + 1))
sd_pap <- df %>% filter(!is.na(mPAP))
sd_pap <- sd(log(sd_pap$mPAP + 1))
sd_bnp <- df %>% filter(!is.na(bnp_level))
sd_bnp <- sd(log(sd_bnp$bnp_level + 1))

#also get sd data for non-log-transformed data
sd_fvc <- df %>% filter(!is.na(FVCp))
sd_fvc <- sd(sd_fvc$FVCp)
sd_egfr <- df %>% filter(!is.na(eGFR))
sd_egfr <- sd(sd_egfr$eGFR)
sd_rap <- df %>% filter(!is.na(mRAP))
sd_rap <- sd(sd_rap$mRAP)

lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'BMI', lm_res_adj$estimate/sd_bmi, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'CRP', lm_res_adj$estimate/sd_crp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'FVCp', lm_res_adj$estimate/sd_fvc, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'mRAP', lm_res_adj$estimate/sd_rap, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'eGFR', lm_res_adj$estimate/sd_egfr, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'PVR', lm_res_adj$estimate/sd_pvr, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'BaseISWD', lm_res_adj$estimate/sd_iswd, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'PAWP', lm_res_adj$estimate/sd_pawp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'mPAP', lm_res_adj$estimate/sd_pap, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'CardiacOutput', lm_res_adj$estimate/sd_co, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'CardiacIndex', lm_res_adj$estimate/sd_ci, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'bnp_level', lm_res_adj$estimate/sd_bnp, lm_res_adj$estimate)

lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'BMI', lm_res_adj$std.error/sd_bmi, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'CRP', lm_res_adj$std.error/sd_crp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'FVCp', lm_res_adj$std.error/sd_fvc, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'mRAP', lm_res_adj$std.error/sd_rap, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'eGFR', lm_res_adj$std.error/sd_egfr, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'PVR', lm_res_adj$std.error/sd_pvr, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'BaseISWD', lm_res_adj$std.error/sd_iswd, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'PAWP', lm_res_adj$std.error/sd_pawp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'mPAP', lm_res_adj$std.error/sd_pap, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'CardiacOutput', lm_res_adj$std.error/sd_co, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'CardiacIndex', lm_res_adj$std.error/sd_ci, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'bnp_level', lm_res_adj$std.error/sd_bnp, lm_res_adj$std.error)

lm_res_adj$std.error <- lm_res_adj$std.error*1.96

#also generate output with adjustment for SD for BMI and CRP
lm_res_adj$std.error <- ifelse(lm_res_adj$term == 'log(bmi$level)', lm_res_adj$std.error*sd_crp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$term == 'log(bmi$BMI)', lm_res_adj$std.error*sd_bmi, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$term == 'log(df2$level)', lm_res_adj$std.error*sd_crp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$term == 'log(df2$BMI)', lm_res_adj$std.error*sd_bmi, lm_res_adj$std.error)

lm_res_adj$estimate <- ifelse(lm_res_adj$term == 'log(bmi$level)', lm_res_adj$estimate*sd_crp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$term == 'log(bmi$BMI)', lm_res_adj$estimate*sd_bmi, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$term == 'log(df2$level)', lm_res_adj$estimate*sd_crp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$term == 'log(df2$BMI)', lm_res_adj$estimate*sd_bmi, lm_res_adj$estimate)

#save for additional plotting
lm_adjusted_output <- lm_res_adj %>% filter(term %in% c('log(df2$level)', 'log(df2$BMI)', 'log(bmi$level)', 'log(bmi$BMI)'))
saveRDS(lm_adjusted_output, 'linear_models_aspire_adjusted.rds')

#=================================================================
#STEP 4
#=================================================================
#run a correlation matrix for outcome variables
#select the cols for the variables in the lm
mat <- df[,c(3, 6, 10:18)]

library(psych)

corr_mat=corr.test(mat, use='pairwise.complete.obs', method="s")
library(corrplot)

corrplot(corr_mat$r, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black", tl.cex = 0.7, tl.offset = 1, number.cex=0.75,
         tl.col = "black", insig = "pch",
         p.mat = corr_mat$p, sig.level = 0.05,
         title = "Spearman correlations of outcome variables")  

#==================================================================
#STEP 5
#==================================================================
bmi_df <- clean_df

pvri <- df_done_filt %>% select(ID, date, PAWP, mPAP, CardiacIndex, CardiacOutput)
pvri$pvri <- (pvri$mPAP - pvri$PAWP)/pvri$CardiacIndex
pvri$tpr <- pvri$mPAP/pvri$CardiacOutput

bmi_df <- merge(pvri[,c(1,5,7,8)], bmi_df, by='ID')

#Make sure BMI does not influence mortality
#make BMI factorial
bmi_df$weight <- ifelse(bmi_df$BMI <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$BMI >=18.5 & bmi_df$BMI <25, 'normal weight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$BMI >=25 & bmi_df$BMI <30, 'pre-obesity', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$BMI >=30 &  bmi_df$BMI <35, 'obesity class-I', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$BMI >=35 &  bmi_df$BMI <40, 'obesity class-II', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$BMI >=40, 'obesity class-III', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)

mdf <- bmi_df %>% filter(!is.na(weight))
mdf$Gender <- as.factor(mdf$Gender)

mdf <- merge(mdf, df[,c(1,22)])
#now generate survival object
sobj2 <- Surv(mdf$YearsDiagDeath, mdf$LifeStatus, type='right')
sfit2 <- survfit(sobj2 ~ weight, data=mdf)

#==================================================================
#cap at 5 years
md2 <- mdf
md2$LifeStatus <- ifelse(md2$YearsDiagDeath >5, 0, md2$LifeStatus)
md2$YearsDiagDeath <- ifelse(md2$YearsDiagDeath >5, 5, md2$YearsDiagDeath)

#get PVR in WU
md2$PVR <- md2$PVR/80
library('cmprsk')

#add in CRP to extra data for cmprisk
md3 <- md2

#now run CoxPH
so1 <- Surv(md2$YearsDiagDeath, md2$LifeStatus, type='right')
so2 <- survfit(so1 ~ weight, data=md2)

#recode weight groups
levels(md2$weight)

new_order <- c("obesity class-III", "obesity class-II", "obesity class-I", 
               "pre-obesity", "normal weight", "underweight")

# Reorder the factor levels
md2$weight <- factor(md2$weight, levels = new_order)

#==================================
#run code to get comparable plot with Priscilla
mdx <- md2

coxmodel5 <- coxph(Surv(YearsDiagDeath, LifeStatus) ~ BMI + AgeAtDiag + Gender + FinalPrimaryPHDiagnosis, data= mdx, x=TRUE)
summary(coxmodel5)

library(contsurvplot)

plot_surv_lines(time="YearsDiagDeath",
                status="LifeStatus",
                variable="BMI",
                data=mdx,
                model=coxmodel5,
                horizon=c(18.5, 25, 30,35,40))


plot_surv_3Dsurface(time="YearsDiagDeath",
                    status="LifeStatus",
                    variable="BMI",
                    data=mdx,
                    model=coxmodel5,
                    interactive = TRUE)

#==================================================

#do KM
p <- ggsurvplot(so2, data=md2, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, xlab='Time since diagnosis (years)', title='survival differences based on BMI', xlim=c(0,5), legend.title='Group', break.x.by=1)
p <- ggsurvplot(so2, data=md2, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, xlab='Time since diagnosis (years)', title='survival differences based on BMI', xlim=c(0,5), legend.title='Group', break.x.by=1, legend.labs=c("obesity class-III", "obesity class-II", "obesity class-I", "pre-obesity", "normal weight", "underweight"))
#also do coxPH

cox_bmi3 <- coxph(so1 ~ BMI + Gender + AgeAtDiag + PVR + mRAP + FinalPrimaryPHDiagnosis, data=md2)
ggforest(cox_bmi3, data=md2)

md2$FinalPrimaryPHDiagnosis <- recode_factor(md2$FinalPrimaryPHDiagnosis, "B - Pulmonary Arterial Hypertension" = "Group 1 PAH", "D - Pulmonary Hypertension due to Left Heart Disease" = "Group 2 PH",
                                             "E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia" = "Group 3 PH",  "F - Chronic Thromboembolic Pulmonary Hypertension" = "Group 4 PH", "G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms" = "Group 5 PH")


cox_bmi3 <- coxph(so1 ~ weight + Gender + AgeAtDiag + PVR + mRAP + FinalPrimaryPHDiagnosis, data=md2)
ggforest(cox_bmi3, data=md2)
summary_cox <- summary(cox_bmi3)
#save this for plotting later
df_coef <- as.data.frame(summary_cox$coefficients)
saveRDS(df_coef, 'Coefficients_ASPIRE_Cox_PH_per_BMI_group_v1_3May24.rds')

#repeat this with scaled vars
md2$pvr_std=scale(md2$PVR)
md2$rap_std=scale(md2$mRAP)
md2$age_std=scale(md2$AgeAtDiag)

#do cox
cox_bmi4 <- coxph(so1 ~ weight + Gender + age_std + pvr_std + rap_std + FinalPrimaryPHDiagnosis, data=md2)
ggforest(cox_bmi4, data=md2)
df_coef2 <- summary(cox_bmi4)
df_coef2 <- as.data.frame(df_coef2$coefficients)
saveRDS(df_coef2, 'SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_v1_3May24.rds')


#get CoxPH per subgroup
results_list <- list()

# Loop through each level of 'FinalPrimaryPHDiagnosis'
for(group in unique(md2$FinalPrimaryPHDiagnosis)) {
  subset_data <- filter(md2, FinalPrimaryPHDiagnosis == group)
  # Run Cox Proportional Hazards model
  cox_model <- coxph(Surv(YearsDiagDeath, LifeStatus) ~ weight + Gender + age_std + pvr_std + rap_std, data = subset_data)
  summary_cox <- summary(cox_model)
  # Extracting coefficients, p-values, and confidence intervals
  cox_res <- summary_cox$coefficients
  results_df <- as.data.frame(cox_res)
  results_list[[group]] <- results_df
  # Print summary (optional)
  print(summary_cox)
} 

# Combine all results into one dataframe
final_results <- do.call(rbind, lapply(names(results_list), function(x) cbind(Group = x, results_list[[x]])))


# Save coefficients data frame
saveRDS(final_results, 'SCALED_Coefficients_All_Groups_PH_per_BMI_v1_7May24.rds')

####################################
#now remove inital worse survival
md_rem <- md2 %>% filter(YearsDiagDeath >=1)
so123 <- Surv(md_rem$YearsDiagDeath, md_rem$LifeStatus, type='right')
so2a <- survfit(so123 ~ weight, data=md_rem)

p <- ggsurvplot(so2a, data=md_rem, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, xlab='Time since diagnosis (years)', title='survival differences based on BMI', xlim=c(0,5), legend.title='Group', break.x.by=1, legend.labs=c("obesity class-III", "obesity class-II", "obesity class-I", "pre-obesity", "normal weight", "underweight"))  

cox_bmi3 <- coxph(so123 ~ weight + Gender + age_std + pvr_std + rap_std + FinalPrimaryPHDiagnosis, data=md_rem)
ggforest(cox_bmi3, data=md_rem)
summary_cox <- summary(cox_bmi3)
df_coef <- as.data.frame(summary_cox$coefficients)
saveRDS(df_coef, 'SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_1yr_removed_v1_7May24.rds')

####################################
#do competing outcomes plot
surv_obj <- with(md3, Surv(time = YearsDiagDeath, event = LifeStatus))

md4 <- md3 %>% select(ID, YearsDiagDeath, LifeStatus, weight, level, FinalPrimaryPHDiagnosis)
md4 <- na.omit(md4)

#now give it another go removing extreme outliers
md4 <- md4 %>% filter(level <=50)

md41 <- md4
md41$weight <- as.character(md41$weight)
md41$weight <- ifelse(md41$weight %in% c('obesity class-I', 'obesity class-II', 'obesity class-III'), 'obese', md41$weight)
md41$weight <- as.factor(md41$weight)

md4a <- md41 %>% filter(FinalPrimaryPHDiagnosis == 'B - Pulmonary Arterial Hypertension')
md4b <- md41 %>% filter(FinalPrimaryPHDiagnosis == 'D - Pulmonary Hypertension due to Left Heart Disease')
md4c <- md41 %>% filter(FinalPrimaryPHDiagnosis == 'E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia')
md4d <- md41 %>% filter(FinalPrimaryPHDiagnosis == 'F - Chronic Thromboembolic Pulmonary Hypertension')
md4e <- md41 %>% filter(FinalPrimaryPHDiagnosis == 'G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms ')

#############################################################
#repeat this with pspline
cox_model_explicit <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4)
c2 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=4) + weight, data = md4)
c3 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=5) + weight, data = md4)
c4 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=6) + weight, data = md4)
c5 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=7) + weight, data = md4)

BIC(cox_model_explicit, c2, c3,c4,c5)

# Check the summary
summary(cox_model_explicit)

m=cox_model_explicit # because I am a lazy typist
termplot(m,term=2,se=TRUE)

levels_weight <- levels(md4$weight)
seq_crp <- seq(from = min(md4$level, na.rm = TRUE), to = max(md4$level, na.rm = TRUE), length.out = 100)
new_data <- data.frame(
  level = rep(seq_crp, times = length(levels_weight)),
  weight = rep(levels_weight, each = length(seq_crp))
)
new_data$YearsDiagDeath <- rep(5, nrow(new_data))  # Assuming 0 is an arbitrary placeholder
new_data$LifeStatus <- rep(0, nrow(new_data))      # Assuming 0 indicates the event has not occurred# Check the structure to ensure it matches expe
new_data$weight <- factor(new_data$weight, levels = levels(md4$weight))# Confirm data structure
str(new_data)
summary(new_data$level)# Try prediction again

fit_surv <- survfit(cox_model_explicit, newdata = new_data)
surv_probs <- fit_surv$surv[length(fit_surv$time),]
new_data$surv_prob <- fit_surv$surv[length(fit_surv$time),]
## surv_times <- fit_surv$time# Extract survival probabilities
## last_indices <- sapply(split(fit_surv$time, f = rep(1:length(fit_surv$strata), each = fit_surv$strata)), max)
## new_data$surv_prob <- surv_probs[last_indices]
new_data$weight <- factor(new_data$weight, levels = c("underweight", "normal weight", "pre-obesity", "obesity class-I", "obesity class-II", "obesity class-III"))

library(ggplot2)
ggplot(new_data, aes(x = level, y = surv_prob, group = weight, color = weight)) +
  geom_line() +  # Adds the lines for each weight group
  labs(
    title = "Survival Probability by diagnostic CRP level and BMI group",
    x = "CRP Level",
    y = "Survival Probability"
  ) +
  theme_minimal() +
  guides(colour = guide_legend(reverse=TRUE)) +
  scale_color_brewer(palette = "Set1")



###############################################
#1 PAH
#repeat per invidual group
#repeat this with pspline
cox_model_explicit1 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4a)
cox_model_explicit2 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=4) + weight, data = md4a)
cox_model_explicit3 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=5) + weight, data = md4a)
cox_model_explicit4 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=6) + weight, data = md4a)
cox_model_explicit5 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=7) + weight, data = md4a)

#compare BIC scores
BIC(cox_model_explicit1)
BIC(cox_model_explicit2)
BIC(cox_model_explicit3)
BIC(cox_model_explicit4)
BIC(cox_model_explicit5)
#3 df is best here (lowest BIC score)
#scores:
#BIC(cox_model_explicit1)
#10508.33
#BIC(cox_model_explicit2)
#10511.35
#BIC(cox_model_explicit3)
#10516.16
#BIC(cox_model_explicit4)
#10520.62
#BIC(cox_model_explicit5)
#10525.2

cox_model_explicit <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4a)
# Check the summary
summary(cox_model_explicit)

m=cox_model_explicit # because I am a lazy typist
termplot(m,term=2,se=TRUE)

levels_weight <- levels(md4a$weight)
seq_crp <- seq(from = min(md4a$level, na.rm = TRUE), to = max(md4a$level, na.rm = TRUE), length.out = 100)
new_data <- data.frame(
  level = rep(seq_crp, times = length(levels_weight)),
  weight = rep(levels_weight, each = length(seq_crp))
)
new_data$YearsDiagDeath <- rep(5, nrow(new_data))  # Assuming 0 is an arbitrary placeholder
new_data$LifeStatus <- rep(0, nrow(new_data))      # Assuming 0 indicates the event has not occurred# Check the structure to ensure it matches expe
new_data$weight <- factor(new_data$weight, levels = levels(md4a$weight))# Confirm data structure
str(new_data)
summary(new_data$level)# Try prediction again

#fit_surv <- survfit(cox_model_explicit, newdata = new_data)
#surv_probs <- fit_surv$surv[length(fit_surv$time),]
#new_data$surv_prob <- fit_surv$surv[length(fit_surv$time),]

pred <- predict(cox_model_explicit, newdata = new_data, type = "lp", se.fit = TRUE)
baseline_hazard <- basehaz(cox_model_explicit, centered = FALSE)

# Calculate the survival probabilities
surv_prob <- exp(-exp(pred$fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])  # Assuming 5 years
new_data$surv_prob <- surv_prob

# Calculate confidence intervals on the linear predictor scale
new_data$lower_ci <- exp(-exp(pred$fit + 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])
new_data$upper_ci <- exp(-exp(pred$fit - 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])


## surv_times <- fit_surv$time# Extract survival probabilities
## last_indices <- sapply(split(fit_surv$time, f = rep(1:length(fit_surv$strata), each = fit_surv$strata)), max)
## new_data$surv_prob <- surv_probs[last_indices]
new_data$weight <- factor(new_data$weight, levels = c("underweight", "normal weight", "pre-obesity", "obese"))

library(ggplot2)
ggplot(new_data, aes(x = level, y = surv_prob, group = weight, color = weight)) +
  geom_line() +  # Adds the lines for each weight group
  labs(
    title = "Survival Probability by diagnostic CRP level and BMI group - group 1 PAH",
    x = "CRP Level",
    y = "Survival Probability"
  ) +
  theme_minimal() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill=weight, colour=NA), alpha = 0.2) +
  guides(colour = guide_legend(reverse=TRUE)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")


#########################################################################
#2PH
cox_model_explicit1 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4b)
cox_model_explicit2 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=4) + weight, data = md4b)
cox_model_explicit3 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=5) + weight, data = md4b)
cox_model_explicit4 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=6) + weight, data = md4b)
cox_model_explicit5 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=7) + weight, data = md4b)

BIC(cox_model_explicit1)
#3964.248
BIC(cox_model_explicit2)
#3968.597
BIC(cox_model_explicit3)
#3969.275
BIC(cox_model_explicit4)
#3977.633
BIC(cox_model_explicit5)
#3981.754
#again best evidence for df=3

cox_model_explicit <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4b)
# Check the summary
summary(cox_model_explicit)

m=cox_model_explicit # because I am a lazy typist
termplot(m,term=2,se=TRUE)

levels_weight <- levels(md4b$weight)
seq_crp <- seq(from = min(md4$level, na.rm = TRUE), to = max(md4b$level, na.rm = TRUE), length.out = 100)
new_data <- data.frame(
  level = rep(seq_crp, times = length(levels_weight)),
  weight = rep(levels_weight, each = length(seq_crp))
)
new_data$YearsDiagDeath <- rep(5, nrow(new_data))  # Assuming 0 is an arbitrary placeholder
new_data$LifeStatus <- rep(0, nrow(new_data))      # Assuming 0 indicates the event has not occurred# Check the structure to ensure it matches expe
new_data$weight <- factor(new_data$weight, levels = levels(md4b$weight))# Confirm data structure
str(new_data)
summary(new_data$level)# Try prediction again

pred <- predict(cox_model_explicit, newdata = new_data, type = "lp", se.fit = TRUE)
baseline_hazard <- basehaz(cox_model_explicit, centered = FALSE)

# Calculate the survival probabilities
surv_prob <- exp(-exp(pred$fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])  # Assuming 5 years
new_data$surv_prob <- surv_prob

# Calculate confidence intervals on the linear predictor scale
new_data$lower_ci <- exp(-exp(pred$fit + 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])
new_data$upper_ci <- exp(-exp(pred$fit - 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])


new_data$weight <- factor(new_data$weight, levels = c("underweight", "normal weight", "pre-obesity", "obese"))

library(ggplot2)
ggplot(new_data, aes(x = level, y = surv_prob, group = weight, color = weight)) +
  geom_line() +  # Adds the lines for each weight group
  labs(
    title = "Survival Probability by diagnostic CRP level and BMI group - group 2 PH",
    x = "CRP Level",
    y = "Survival Probability"
  ) +
  theme_minimal() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill=weight, colour=NA), alpha = 0.2) +
  guides(colour = guide_legend(reverse=TRUE)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

#############################################################################
#group 3
#repeat this with pspline
cox_model_explicit1 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4c)
cox_model_explicit2 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=4) + weight, data = md4c)
cox_model_explicit3 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=5) + weight, data = md4c)
cox_model_explicit4 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=6) + weight, data = md4c)
cox_model_explicit5 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=7) + weight, data = md4c)

BIC(cox_model_explicit1)
#4962.15
BIC(cox_model_explicit2)
#4965.018
BIC(cox_model_explicit3)
#4969.051
BIC(cox_model_explicit4)
#4973.187
BIC(cox_model_explicit5)
#4977.1

#again best BIC at 3

cox_model_explicit <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4c)
# Check the summary
summary(cox_model_explicit)

m=cox_model_explicit # because I am a lazy typist
termplot(m,term=2,se=TRUE)

levels_weight <- levels(md4c$weight)
seq_crp <- seq(from = min(md4c$level, na.rm = TRUE), to = max(md4c$level, na.rm = TRUE), length.out = 100)
new_data <- data.frame(
  level = rep(seq_crp, times = length(levels_weight)),
  weight = rep(levels_weight, each = length(seq_crp))
)
new_data$YearsDiagDeath <- rep(5, nrow(new_data))  # Assuming 0 is an arbitrary placeholder
new_data$LifeStatus <- rep(0, nrow(new_data))      # Assuming 0 indicates the event has not occurred# Check the structure to ensure it matches expe
new_data$weight <- factor(new_data$weight, levels = levels(md4c$weight))# Confirm data structure
str(new_data)
summary(new_data$level)# Try prediction again

pred <- predict(cox_model_explicit, newdata = new_data, type = "lp", se.fit = TRUE)
baseline_hazard <- basehaz(cox_model_explicit, centered = FALSE)

# Calculate the survival probabilities
surv_prob <- exp(-exp(pred$fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])  # Assuming 5 years
new_data$surv_prob <- surv_prob

# Calculate confidence intervals on the linear predictor scale
new_data$lower_ci <- exp(-exp(pred$fit + 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])
new_data$upper_ci <- exp(-exp(pred$fit - 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])

new_data$weight <- factor(new_data$weight, levels = c("underweight", "normal weight", "pre-obesity", "obese"))

library(ggplot2)
ggplot(new_data, aes(x = level, y = surv_prob, group = weight, color = weight)) +
  geom_line() +  # Adds the lines for each weight group
  labs(
    title = "Survival Probability by diagnostic CRP level and BMI group - group 3 PH",
    x = "CRP Level",
    y = "Survival Probability"
  ) +
  theme_minimal() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill=weight, colour=NA), alpha = 0.2) +
  guides(colour = guide_legend(reverse=TRUE)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")


########################################################################3
#group 4PH
#repeat this with pspline
cox_model_explicit1 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4d)
cox_model_explicit2 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=4) + weight, data = md4d)
cox_model_explicit3 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=5) + weight, data = md4d)
cox_model_explicit4 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=6) + weight, data = md4d)
cox_model_explicit5 <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=7) + weight, data = md4d)

BIC(cox_model_explicit1)
#2489.684
BIC(cox_model_explicit2)
#2493.407
BIC(cox_model_explicit3)
#2497.272
BIC(cox_model_explicit4)
#2501.152
BIC(cox_model_explicit5)
#2505.402

cox_model_explicit <- coxph(Surv(time = YearsDiagDeath, event = LifeStatus) ~ pspline(level, df=3) + weight, data = md4d)

# Check the summary
summary(cox_model_explicit)

m=cox_model_explicit # because I am a lazy typist
termplot(m,term=2,se=TRUE)

levels_weight <- levels(md4d$weight)
seq_crp <- seq(from = min(md4d$level, na.rm = TRUE), to = max(md4d$level, na.rm = TRUE), length.out = 100)
new_data <- data.frame(
  level = rep(seq_crp, times = length(levels_weight)),
  weight = rep(levels_weight, each = length(seq_crp))
)
new_data$YearsDiagDeath <- rep(5, nrow(new_data))  # Assuming 0 is an arbitrary placeholder
new_data$LifeStatus <- rep(0, nrow(new_data))      # Assuming 0 indicates the event has not occurred# Check the structure to ensure it matches expe
new_data$weight <- factor(new_data$weight, levels = levels(md4d$weight))# Confirm data structure
str(new_data)
summary(new_data$level)# Try prediction again

pred <- predict(cox_model_explicit, newdata = new_data, type = "lp", se.fit = TRUE)
baseline_hazard <- basehaz(cox_model_explicit, centered = FALSE)

# Calculate the survival probabilities
surv_prob <- exp(-exp(pred$fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])  # Assuming 5 years

new_data$weight <- factor(new_data$weight, levels = c("underweight", "normal weight", "pre-obesity", "obese"))

# Calculate confidence intervals on the linear predictor scale
new_data$lower_ci <- exp(-exp(pred$fit + 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])
new_data$upper_ci <- exp(-exp(pred$fit - 1.96 * pred$se.fit) * baseline_hazard$hazard[match(5, baseline_hazard$time)])

new_data$weight <- factor(new_data$weight, levels = c("underweight", "normal weight", "pre-obesity", "obesity class-I", "obesity class-II", "obesity class-III"))

library(ggplot2)
ggplot(new_data, aes(x = level, y = surv_prob, group = weight, color = weight)) +
  geom_line() +  # Adds the lines for each weight group
  labs(
    title = "Survival Probability by diagnostic CRP level and BMI group - group 4 PH",
    x = "CRP Level",
    y = "Survival Probability"
  ) +
  theme_minimal() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill=weight, colour=NA), alpha = 0.2) +
  guides(colour = guide_legend(reverse=TRUE)) +
  scale_color_brewer(palette = "Set1")


########################################################################
#group 5
#not enough data to model
#===============================================================

#=================================================================

#now run cox-ph for age + sex + bmi
#relevel first
levels(mdf$FinalPrimaryPHDiagnosis)
mdf$FinalPrimaryPHDiagnosis <- recode_factor(mdf$FinalPrimaryPHDiagnosis, "B - Pulmonary Arterial Hypertension" = "Group 1 PAH", "D - Pulmonary Hypertension due to Left Heart Disease" = "Group 2 PH",
                                             "E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia" = "Group 3 PH",  "F - Chronic Thromboembolic Pulmonary Hypertension" = "CTEPH", "G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms" = "Group 5 PH")


library(Publish)
uv1 <- univariateTable(weight ~ Gender + AgeAtDiag + FinalPrimaryPHDiagnosis + BMI + level + PVR + pvri + tpr + CardiacIndex + CardiacOutput + mRAP + mPAP + PAWP + FVCp + eGFR + BaseISWD, data=mdf, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE, summary.format = 'median(x) [iqr(x)]')
uv2 <- summary(uv1)

uv3 <- uv2[,c(1,2,8,3,7,4:6, 9, 10)]

hist(mdf$AgeAtDiag)
hist(mdf$BMI)
hist(mdf$level)
hist(mdf$pvri)
hist(mdf$tpr)
hist(mdf$CardiacIndex)
hist(mdf$CardiacOutput)
hist(mdf$mPAP)
hist(mdf$mRAP)
hist(mdf$PAWP)
hist(mdf$FVCp)
hist(mdf$eGFR)
hist(mdf$BaseISWD)

library(broom)

do_kruskal = function(myvar) {
 output_krusk <- kruskal.test(mdf[[myvar]] ~ mdf$weight)
 output_krusk = tidy(output_krusk)
 output_krusk$var <- myvar
 return(output_krusk)
}

do_kruskal('AgeAtDiag')
#this works
myvar = c('AgeAtDiag', 'BMI', 'level', 'PVR', 'pvri', 'tpr',  'CardiacIndex', 'CardiacOutput', 'mRAP', 'mPAP', 'PAWP', 'FVCp', 'eGFR', 'BaseISWD')
o1 <- lapply(myvar, do_kruskal)
o2 <- bind_rows(o1)
uv3$num <- 1:nrow(uv3)

uv4 <- merge(uv3, o2[,c(5,2)], by.x='Variable', by.y='var', all=T)
uv4 <- uv4[order(uv4$num),]

write.csv2(uv4, file='obesity_baseline_differences.csv')
#===========================================================
#STEP 4
#===========================================================
#Check survival differences between CRP groups
#take DIAGNOSTIC or FIRST RECORDED CRP
clean_df$crp_above_5 <- ifelse(clean_df$level >5, 'yes', 'no')


#######################################################################
#be stringent, set IQR of 4 as threshold for later filtering
iqr_4 <- 4*IQR(clean_df$level) + median(clean_df$level)
iqr_4
#IQR of 4 corresponds to a CRP of 40.8
excl_crp <- clean_df %>% filter(level <iqr_4)
#########################################################################

clean_df$crp_above_5 <- as.factor(clean_df$crp_above_5)

sobj <- Surv(clean_df$YearsDiagDeath, clean_df$LifeStatus, type='right')
sfit <- survfit(sobj ~ crp_above_5, data=clean_df)

#=================================================================
#run with capping survival at 5 years
cdf <- clean_df
cdf$LifeStatus <- ifelse(cdf$YearsDiagDeath >5, 0, cdf$LifeStatus)
cdf$YearsDiagDeath <- ifelse(cdf$YearsDiagDeath >5, 5, cdf$YearsDiagDeath)

#now run CoxPH
soc1 <- Surv(cdf$YearsDiagDeath, cdf$LifeStatus, type='right')
soc2 <- survfit(soc1 ~ crp_above_5, data=cdf)

#do KM
p <- ggsurvplot(soc2, data=cdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, xlab='Time since diagnosis (years)', title='survival differences based on diagnostic CRP', xlim=c(0,5), legend.title='CRP level', break.x.by=1, legend.labs=c('CRP >5', 'CRP <=5'), palette = c('darkred', 'darkblue'))
#also do coxPH

#get PVR in WU
cdf$PVR <- cdf$PVR/80

cox5 <- coxph(soc1 ~ crp_above_5 +Gender + BMI + AgeAtDiag + mRAP + PVR + FinalPrimaryPHDiagnosis, data=cdf)
ggforest(cox5, data=cdf)
sdf <- summary(cox5)
sdf2 <- as.data.frame(sdf$coefficients)
saveRDS(sdf2, 'coefficient_cox_ph_CRP_groups_aspire_v1_7May.rds')

cdf$FinalPrimaryPHDiagnosis <- recode_factor(cdf$FinalPrimaryPHDiagnosis, "B - Pulmonary Arterial Hypertension" = "Group 1 PAH", "D - Pulmonary Hypertension due to Left Heart Disease" = "Group 2 PH",
                                             "E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia" = "Group 3 PH",  "F - Chronic Thromboembolic Pulmonary Hypertension" = "CTEPH", "G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms" = "Group 5 PH")


#repeat this for scaled vars
cdf$sc_pvr <- scale(cdf$PVR)
cdf$sc_bmi <- scale(cdf$BMI)
cdf$sc_rap <- scale(cdf$mRAP)
cdf$sc_age <- scale(cdf$AgeAtDiag)

cdf2 <- cdf

cox5 <- coxph(soc1 ~ crp_above_5 +Gender + sc_bmi + sc_age + sc_rap + sc_pvr + FinalPrimaryPHDiagnosis, data=cdf)
ggforest(cox5, data=cdf)

s5 <- summary(cox5)
s51 <- as.data.frame(s5$coefficients)

saveRDS(s51, file='SCALED_coefficients_coxPHallPH_per_CRP.rds')

#do KM
p <- ggsurvplot(soc2, data=cdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, xlab='Time since diagnosis (years)', title='survival differences based on diagnostic CRP', xlim=c(0,5), legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'))


#repeat capped survival for >4IQR removed data

cdf <- excl_crp
cdf$LifeStatus <- ifelse(cdf$YearsDiagDeath >5, 0, cdf$LifeStatus)
cdf$YearsDiagDeath <- ifelse(cdf$YearsDiagDeath >5, 5, cdf$YearsDiagDeath)

#now run CoxPH
soc1 <- Surv(cdf$YearsDiagDeath, cdf$LifeStatus, type='right')
soc2 <- survfit(soc1 ~ crp_above_5, data=cdf)

#do KM
p <- ggsurvplot(soc2, data=cdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, xlab='Time since diagnosis (years)', title='survival differences based on diagnostic CRP, outliers removed', xlim=c(0,5), legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'))

#=================================================================
#do KM per disease group
#and do CoxPH all censored at 5 years post diagnosis

cdf <- cdf2

#cdf$FinalPrimaryPHDiagnosis <- recode_factor(cdf$FinalPrimaryPHDiagnosis, "B - Pulmonary Arterial Hypertension" = "Group 1 PAH", "D - Pulmonary Hypertension due to Left Heart Disease" = "Group 2 PH",
                                                  "E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia" = "Group 3 PH",  "F - Chronic Thromboembolic Pulmonary Hypertension" = "CTEPH", "G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms" = "Group 5 PH")


d1 <- cdf %>% filter(FinalPrimaryPHDiagnosis == 'Group 1 PAH')
d2 <- cdf %>% filter(FinalPrimaryPHDiagnosis == 'Group 2 PH')
d3 <- cdf %>% filter(FinalPrimaryPHDiagnosis == 'Group 3 PH')
d4 <- cdf %>% filter(FinalPrimaryPHDiagnosis == 'CTEPH')
d5 <- cdf %>% filter(FinalPrimaryPHDiagnosis == 'Group 5 PH') 

sa <- Surv(d1$YearsDiagDeath, d1$LifeStatus, type='right')
s2 <- survfit(sa ~ crp_above_5, data=d1)

sb <- Surv(d2$YearsDiagDeath, d2$LifeStatus, type='right')
s3 <- survfit(sb ~ crp_above_5, data=d2)

sc <- Surv(d3$YearsDiagDeath, d3$LifeStatus, type='right')
s4 <- survfit(sc ~ crp_above_5, data=d3)

sd <- Surv(d4$YearsDiagDeath, d4$LifeStatus, type='right')
s5 <- survfit(sd ~ crp_above_5, data=d4)

se <- Surv(d5$YearsDiagDeath, d5$LifeStatus, type='right')
s6 <- survfit(se ~ crp_above_5, data=d5)

pa <- ggsurvplot(s2, data=d1, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences for CRP in PAH', xlim=c(0,5), xlab= 'Time (years since diagnosis)', legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'), theme_risktable = theme(axis.text = element_text(size = 25)))
pb <- ggsurvplot(s3, data=d2, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences for CRP in Group 2 PH', xlim=c(0,5), xlab= 'Time (years since diagnosis)', legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'), theme_risktable = theme(axis.text = element_text(size = 25)))
pc<- ggsurvplot(s4, data=d3, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences for CRP in Group 3 PH', xlim=c(0,5), xlab= 'Time (years since diagnosis)', legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'), theme_risktable = theme(axis.text = element_text(size = 25)))
pd <- ggsurvplot(s5, data=d4, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences for CRP in Group 4 PH', xlim=c(0,5), xlab= 'Time (years since diagnosis)', legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'), theme_risktable = theme(axis.text = element_text(size = 25)))
pe <- ggsurvplot(s6, data=d5, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences for CRP in Group 5 PH', xlim=c(0,5), xlab= 'Time (years since diagnosis)', legend.title='CRP level', break.x.by=1, legend.labs=c('CRP <=5', 'CRP >5'), palette = c('darkblue', 'darkred'), theme_risktable = theme(axis.text = element_text(size = 25)))


combplot <- arrange_ggsurvplots(list(pa, pb, pc, pd, pe), ncol=2, nrow = 3)
ggsave('combined_plot.pdf', combplot, widt=15, height=15)

#do coxmodels
cox1 <- coxph(sa ~ crp_above_5 +Gender + sc_bmi + sc_age + sc_rap + sc_pvr, data=d1)
cox2 <- coxph(sb ~ crp_above_5 +Gender + sc_bmi + sc_age + sc_rap + sc_pvr, data=d2)
cox3 <- coxph(sc ~ crp_above_5 +Gender + sc_bmi + sc_age + sc_rap + sc_pvr, data=d3)
cox4 <- coxph(sd ~ crp_above_5 +Gender + sc_bmi + sc_age + sc_rap + sc_pvr, data=d4)
cox5 <- coxph(se ~ crp_above_5 +Gender + sc_bmi + sc_age + sc_rap + sc_pvr, data=d5)

sum1 <- summary(cox1)
sum1 <- as.data.frame(sum1$coefficients)
sum1$group <- 'Group 1 PAH'
#saveRDS(sum1, 'scaled_summary_PAH_cox_aspire.rds')

sum2 <- summary(cox2)
sum2 <- as.data.frame(sum2$coefficients)
sum2$group <- 'Group 2 PH'
#saveRDS(sum2, 'scaled_summary_PH2_cox_aspire.rds')

sum3 <- summary(cox3)
sum3 <- as.data.frame(sum3$coefficients)
sum3$group <- 'Group 3 PH'
#saveRDS(sum3, 'scaled_summary_PH3_cox_aspire.rds')

sum4 <- summary(cox4)
sum4 <- as.data.frame(sum4$coefficients)
sum4$group <- 'Group 4 PAH'
#saveRDS(sum4, 'scaled_summary_PH4_cox_aspire.rds')

sum5 <- summary(cox5)
sum5 <- as.data.frame(sum5$coefficients)
sum5$group <- 'Group 5 PH'
#saveRDS(sum5, 'scaled_summary_PH5_cox_aspire.rds')

#merge all
coefs <- rbind(sum1, sum2, sum3, sum4, sum5)

saveRDS(coefs, file='coefficients_SCALED_effects_cox_crp_per_PH_subgroup.rds')
#================================================
#STEP 5
#=================================================
#now compare clinical differences of interest between groups
clin123 = clean_df

#first plot clinical differences
labs <- c('<=5', '>5')
box_fills <- c('darkblue', 'darkred')

library(rstatix)

a1 <- ggplot(clin123, aes(x=crp_above_5, y=PAWP)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('PAWP') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('Wedge pressure (mmHg)') +
  xlab('CRP group') +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)


a2 <- ggplot(clin123, aes(x=crp_above_5, y=mRAP)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('RAP') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('RAP (mmHg)') +
  xlab('CRP group') +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a3 <- ggplot(clin123, aes(x=crp_above_5, y=BMI)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('BMI') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('BMI (kg/m^2)') +
  xlab('CRP group') +
 # stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a4 <- ggplot(clin123, aes(x=crp_above_5, y=eGFR)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('eGFR') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('eGFR (ml/min/1.73m^2)') +
  xlab('CRP group') +
#  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)


a5 <- ggplot(clin123, aes(x=crp_above_5, y=BaseISWD)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Shuttle walk') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD (metres)') +
  xlab('CRP group') +
 # stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a6 <- ggplot(clin123, aes(x=crp_above_5, y=FVCp)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('FVC') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('FVC (%pred)') +
  xlab('CRP group') +
  #stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a7 <- ggplot(clin123, aes(x=crp_above_5, y=PVR)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('PVR') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('PVR (Dynes/sec)') +
  xlab('CRP group') +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs) +
  theme_bw() + 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)
  
#combine all
library(patchwork)
  
a1 + a2 + a3 + a4 +a5 + a6 + a7 + plot_layout(ncol=7)

#######################################################
#now also per subtype

clin123$FinalPrimaryPHDiagnosis <- recode_factor(clin123$FinalPrimaryPHDiagnosis, "B - Pulmonary Arterial Hypertension" = "Group 1 PAH", "D - Pulmonary Hypertension due to Left Heart Disease" = "Group 2 PH",
                                             "E - Pulmonary Hypertension due to Lung Disease and/or Hypoxia" = "Group 3 PH",  "F - Chronic Thromboembolic Pulmonary Hypertension" = "CTEPH", "G - Pulmonary Hypertension with Unclear/Multifactorial Mechanisms" = "Group 5 PH")


p1<-  ggplot(clin123, aes(x=crp_above_5, y=PAWP)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('PAWP') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('Wedge pressure (mmHg)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))
  
  
p2 <-  ggplot(clin123, aes(x=crp_above_5, y=mRAP)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('RAP') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('RAP (mmHg)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))
  
p3 <-   ggplot(clin123, aes(x=crp_above_5, y=BMI)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('BMI') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('BMI (kg/m^2)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))
  
p4 <-  ggplot(clin123, aes(x=crp_above_5, y=eGFR)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('eGFR') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('eGFR (ml/min/1.73m^2)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))
  
  
p5 <-  ggplot(clin123, aes(x=crp_above_5, y=BaseISWD)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('Shuttle walk') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD (metres)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))
  
p6 <-  ggplot(clin123, aes(x=crp_above_5, y=FVCp)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('FVC') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('FVC (%pred)') +
  xlab('CRP group') +
  stat_compare_means(method = "t.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))

p7 <-  ggplot(clin123, aes(x=crp_above_5, y=PVR)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=2) +
  ggtitle('PVR') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('PVR (Dynes/sec)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test", vjust = 1, label = 'p.format') +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)+ 
    scale_y_continuous(expand=expansion(mult = c(0, .1)))

p1 + p2 + p3 + p4 +p5 + p6 + p7 + plot_layout(ncol=4)

#####################################################################3

clin123$who_functional_class <- as.factor(clin123$who_functional_class)

chisq.test(clin123$who_functional_class, clin123$crp_above_5)

#get CI
pvri <- df_done_filt %>% select(ID, date, PAWP, mPAP, CardiacIndex, CardiacOutput)
pvri$pvri <- (pvri$mPAP - pvri$PAWP)/pvri$CardiacIndex
pvri$tpr <- pvri$mPAP/pvri$CardiacOutput

clin1234 <- merge(clin123[,c(1,7,21)], pvri, by='ID')
hist(clin1234$pvri)

abc <- merge(bmi_df, pvri, by='ID')

ggplot(clin1234, aes(x=crp_above_5, y=pvri)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in PVR-I between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('Indexed PVR') +
  xlab('CRP group') +
#  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw() +
  facet_wrap(~FinalPrimaryPHDiagnosis)

wilcox.test(clin123$BMI ~ clin123$crp_above_5)
wilcox.test(clin123$PAWP ~ clin123$crp_above_5)
wilcox.test(clin123$mRAP ~ clin123$crp_above_5)
wilcox.test(clin123$FVCp ~ clin123$crp_above_5)
wilcox.test(clin123$eGFR ~ clin123$crp_above_5)
wilcox.test(clin123$PVR ~ clin123$crp_above_5)

crp2 <- merge(df_done3, crp_df, by.x='sth_id', by.y='ID')
crp2$who_functional_class <- as.factor(crp2$who_functional_class)
chisq.test(crp2$who_functional_class, crp2$above_5)
wilcox.test(crp2$walking_distance ~ crp2$above_5)

#we can discuss adding in extra vars and pval corrections

#plot CRP trajectory per group! 

clean_df

plot_df <- df_crp
plot_df <- merge(plot_df, clean_df[,c(1,21)], by='ID')
plot_df$time_since_diagnosis <- as.numeric(plot_df$date - plot_df$DateofFinalDiagnosis)/365.242199
#remove values <-1 years pre-diagnosis
plot_df <- plot_df %>% filter(time_since_diagnosis >=-1)

#round so that we can include half years
library(plyr)
plot_df$time_since_diagnosis <- round_any(plot_df$time_since_diagnosis, 0.5)
#also cut at 10 years
plot_df <- plot_df %>% filter(time_since_diagnosis <=10)

#do some testing if trajectory first 5 years differs significantly
plot2 <- plot_df %>% filter(time_since_diagnosis <=5)
plot2 <- plot2[,c(1,3,21,22)]
crp1 <- plot2
crp1$time_since_diagnosis <- as.numeric(crp1$time_since_diagnosis)
#remove diagnostic visit
crp1 <- crp1 %>% group_by(ID) %>% filter(time_since_diagnosis !=0)
crp1$ID <- as.factor(crp1$ID)
crp1$crp_above_5 <- as.factor(crp1$crp_above_5)
crp1$time_since_diagnosis <- as.numeric(crp1$time_since_diagnosis)

#get median value per half year
crp2 <- crp1 %>% dplyr::group_by(ID, time_since_diagnosis, crp_above_5) %>% dplyr::summarise(median_crp_time = median(level))

#get 2 way repeated measure time series
#first log-transform
crp2$log_med_crp <- log(crp2$median_crp_time)

#this doesn't work --> try something else:
model.aov <- aov(crp2$median_crp_time ~ 
                   crp2$crp_above_5 * crp2$time_since_diagnosis + 
                   Error(crp2$ID/(crp2$crp_above_5*crp2$time_since_diagnosis)))

summary(model.aov)
#======================================================

basic_plot <- plot_df %>% dplyr::group_by(crp_above_5, time_since_diagnosis) %>% dplyr::summarise(mean = mean(level), stdev = sd(level), meds = median(level), iqr = 0.5*IQR(level))

ggplot(data=basic_plot, aes(x=time_since_diagnosis, y=meds, colour=as.factor(crp_above_5))) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=meds-iqr, ymax=meds+iqr), width=.2) +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5,6,7,8,9,10)) +
  ggtitle('CRP levels in all cause PH over time (median + IQR)') +
  ylab('Median CRP level (mg/l)') +
  xlab('Time since diagnosis (years)') +
  theme_bw() +
  scale_color_manual(values=c("darkblue", "darkred")) +
  guides(color = guide_legend(title = "Diagnostic CRP level >5"))


#get geometric means
library(plotrix)
gm_crp <- plot_df %>% dplyr::group_by(crp_above_5, time_since_diagnosis) %>% dplyr::summarise(geo_mean = exp(mean(log(level))), ci = plotrix::std.error(level))

plot_gm <- ggplot(gm_crp, aes(x=time_since_diagnosis, y=geo_mean, colour=as.factor(crp_above_5))) +
  geom_point() +
  geom_line() +
  xlab("Time since diagnosis (years)") +
  ylab("CRP level (mg/ml)") +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5,6,7,8,9,10)) +
  geom_errorbar(aes(ymin=geo_mean-ci, ymax=geo_mean+ci), width=.2) +
  ggtitle("Geometric means for CRP, with 95% CI per group") + 
  scale_color_manual(values=c("darkblue", "darkred")) +
  theme_bw() +
  guides(color = guide_legend(title = "Diagnostic CRP level >5"))

#cut at 5 years!
gm_crp2 <- gm_crp %>% filter(time_since_diagnosis <=5)

plot_gm <- ggplot(gm_crp2, aes(x=time_since_diagnosis, y=geo_mean, colour=as.factor(crp_above_5))) +
  geom_point() +
  geom_line() +
  xlab("Time since diagnosis (years)") +
  ylab("CRP level (mg/ml)") +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5,6,7,8,9,10)) +
  geom_errorbar(aes(ymin=geo_mean-ci, ymax=geo_mean+ci), width=.2) +
  ggtitle("Geometric means for CRP, with 95% CI per group - ASPIRE") + 
  scale_color_manual(values=c("darkblue", "darkred")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 20)) +
  guides(color = guide_legend(title = "CRP >=5"))


#check baseline differences
library(Publish)
clean_df$FinalPrimaryPHDiagnosis <- as.factor(clean_df$FinalPrimaryPHDiagnosis)
clean_df$who_functional_class <- as.factor(clean_df$who_functional_class)
a1 <-univariateTable(FinalPrimaryPHDiagnosis ~ Gender + AgeAtDiag + BMI + level + mRAP + mPAP + PAWP + CardiacOutput + PVR + BaseISWD + who_functional_class, data=clean_df, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
a2 <- summary(a1)

#check other baseline differences between groups
b1 <-univariateTable(crp_above_5 ~ Gender + AgeAtDiag + BMI + level + mRAP + mPAP + PAWP + CardiacOutput + PVR + BaseISWD + who_functional_class + FinalPrimaryPHDiagnosis, data=clean_df, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
b2 <- summary(b1)


#======================================================================
#EdB continued here on 22-03-2024
#assess treatment effects
iswd <- df_done_filt %>% select(ID, BMI, BaseISWD, FUISWD, BaselineE10score, FUE10atanytimeScore, FinalPrimaryPHDiagnosis, AgeatBaseISWD, AgeatFUISWT)

#select other groups too
crp_group <- clean_df %>% select(ID, crp_above_5)
bmi_group <- bmi_df %>% select(ID, weight)

groups <- merge(crp_group, bmi_group, by='ID')
groups$crp_above_5 <- as.factor(groups$crp_above_5)

iswd <- merge(groups, iswd, by='ID')
#calculate ISWD difference
iswd$difference <- iswd$FUISWD-iswd$BaseISWD 
iswd$percent_diff <- (iswd$difference/iswd$BaseISWD)*100

#plot distributions
hist(iswd$difference)
#roughly normal
hist(iswd$percent_diff)
#percentual difference is inf when started with 0
iswd$percent_diff <- ifelse(iswd$percent_diff == Inf, iswd$difference, iswd$percent_diff)
#same when this difference is 0
iswd$percent_diff <- ifelse(iswd$difference == 0, 0, iswd$percent_diff)
hist(iswd$percent_diff)
hist(scale(log(iswd$percent_diff + 101)))


#regress out baseline ISWD to see effect
hist(log(iswd$BaseISWD +5))
#log transformation is required to make baseISWD normal
regression_test <- iswd %>% filter(!is.na(difference))
regression_test$bl_log_iswd <- log(regression_test$BaseISWD +5)
#now regress out this
fit_bl <- lm(regression_test$difference ~ regression_test$bl_log_iswd)
regression_test$resids <- fit_bl$residuals
hist(regression_test$resids)
cor.test(regression_test$bl_log_iswd, regression_test$difference)

#remove people without weight
regression_test <- regression_test %>% filter(!is.na(weight))

ggplot(regression_test, aes(x=weight, y=resids)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in residuals of ISWD improvement between weight groups') +
  ylab('ISWD improvement') +
  theme_bw() +
  xlab('BMI group') +
  stat_compare_means(method = "anova")

  
#also propensity match on baseline ISWD
r1 <- regression_test %>% filter(weight == 'normal weight')
r1$class <- 'normal_weight'
r2 <- regression_test %>% filter(weight %in% c('obesity class-I', 'obesity class-II', 'obesity class-III'))
r2$class <- 'obese'

mr <- rbind(r1, r2)
mr$class <- as.factor(mr$class)
levels(mr$class)

library(MatchIt)
matchset <- matchit(formula = class ~ BaseISWD, data= mr, method = 'optimal', distance = 'glm', replace = FALSE, ratio = 1)
mat2 <- MatchIt::match.data(matchset)

#remove horrible matches
mat1 <- mat2 %>% filter(distance <0.55)
#only retaining 136 pts
ggplot(mat1, aes(x=weight, y=difference)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in ISWD improvement between matched BMI groups') +
  ylab('ISWD improvement (m)') +
  xlab('BMI group') +
  stat_compare_means(method = "anova")


#now regress out BMI of difference and %diff
iswd <- iswd %>% filter(!is.na(difference))
fit <- lm(iswd$difference ~ log(iswd$BMI))
bmi_iswd <- iswd %>% filter(!is.na(BMI))
bmi_iswd$residuals_fit <- fit$residuals
hist(bmi_iswd$residuals_fit)
#roughly normal again

#repeat this for percentual differce
fit2 <- lm(bmi_iswd$log_per_diff ~ log(bmi_iswd$BMI))
bmi_iswd$residuals_fit2 <- fit2$residuals
hist(bmi_iswd$residuals_fit2)
#also not really normal...

#plot differences
ggplot(bmi_iswd, aes(x=crp_above_5, y=residuals_fit)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in residuals of ISWD improvement between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD improvement') +
  xlab('CRP group') +
  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs)

ggplot(bmi_iswd, aes(x=crp_above_5, y=residuals_fit2)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in residuals of % ISWD improvement between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD improvement') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs)

#now assess normal differences
ggplot(iswd, aes(x=crp_above_5, y=difference)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in ISWD improvement between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD improvement (m)') +
  xlab('CRP group') +
  theme_bw() +
  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs)

ggplot(iswd, aes(x=crp_above_5, y=percent_diff)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in % ISWD improvement between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD improvement (%)') +
  xlab('CRP group') +
  stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs)

#repeat for BMI groups
levels(bmi_iswd$weight)

#relevel
desired_order <- c("underweight", "normal weight", "pre-obesity", "obesity class-I", "obesity class-II", "obesity class-III")

# Reorder the levels of the factor
bmi_iswd$weight <- factor(bmi_iswd$weight, levels = desired_order)

# Check the new levels to confirm the order
levels(bmi_iswd$weight)

ggplot(bmi_iswd, aes(x=weight, y=difference)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in ISWD improvement between BMI groups in all cause PH') +
  ylab('ISWD improvement (m)') +
  xlab('BMI group') +
  theme_bw() +
  stat_compare_means(method = "anova")

b2 <- bmi_iswd %>% filter(FinalPrimaryPHDiagnosis == 'B - Pulmonary Arterial Hypertension')
b2$agediff <- b2$AgeatFUISWT - b2$AgeatBaseISWD

ggplot(b2, aes(x=weight, y=difference)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in ISWD improvement between BMI groups for group 1 PAH') +
  ylab('ISWD improvement (m)') +
  xlab('BMI group') +
  theme_bw() +
  stat_compare_means(method = "anova")

b3 <- b2 %>% filter(difference >=0)

ggplot(b3, aes(x=weight, y=difference)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in ISWD improvement between BMI groups for group 1 PH') +
  ylab('ISWD improvement (m)') +
  xlab('BMI group') +
  stat_compare_means(method = "anova")

ggplot(bmi_iswd, aes(x=weight, y=percent_diff)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in % ISWD improvement between BMI groups') +
  ylab('ISWD improvement (%)') +
  #ylim(c(-101, 1000)) +
  xlab('BMI group') +
  stat_compare_means(method = "kruskal.test")

#check score differences
iswd$score_diff <- iswd$FUE10atanytimeScore - iswd$BaselineE10score
hist(iswd$score_diff)
iswd$scorediff_perc <- (iswd$score_diff/iswd$BaselineE10score)*100
hist(log(iswd$scorediff_perc + 10))
#log-transforming helps make this parametric

sdf <- iswd %>% filter(!is.na(score_diff))
sdf$logsdf <- log(sdf$scorediff_perc + 101)
hist(sdf$logsdf)

sdf_bmi <- sdf %>% filter(!is.na(weight))


#repeat for BMI groups
ggplot(sdf_bmi, aes(x=weight, y=score_diff)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in E10 score improvement between BMI groups') +
  ylab('ISWD improvement (m)') +
  xlab('BMI group') +
  stat_compare_means(method = "anova")

ggplot(sdf_bmi, aes(x=weight, y=logsdf)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in log % E10 score improvement between BMI groups') +
  ylab('ISWD improvement (%)') +
  #ylim(c(-101, 1000)) +
  xlab('BMI group') +
  stat_compare_means(method = "anova")


#now assess normal differences
ggplot(sdf, aes(x=crp_above_5, y=score_diff)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in ISWD improvement between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD improvement (m)') +
  xlab('CRP group') +
  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs)

ggplot(sdf, aes(x=crp_above_5, y=logsdf)) + 
  geom_boxplot(aes(fill = crp_above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in % ISWD improvement between CRP groups') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('ISWD improvement (%)') +
  xlab('CRP group') +
  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs)



#=====================================================================
#make table for CRP groups
av1 <- univariateTable(crp_above_5 ~ Gender + AgeAtDiag + FinalPrimaryPHDiagnosis + BMI + level + PVR + CardiacOutput + mRAP + mPAP + PAWP + FVCp + eGFR + BaseISWD, data=clin123, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE, summary.format = 'median(x) [iqr(x)]')
av2 <- summary(av1)

library(broom)

do_kruskal2 = function(myvar) {
  output_krusk <- kruskal.test(clin123[[myvar]] ~ clin123$crp_above_5)
  output_krusk = tidy(output_krusk)
  output_krusk$var <- myvar
  return(output_krusk)
}

do_kruskal2('AgeAtDiag')

#this works
myvar = c('AgeAtDiag', 'BMI', 'level', 'PVR', 'CardiacOutput', 'mRAP', 'mPAP', 'PAWP', 'FVCp', 'eGFR', 'BaseISWD')
o1 <- lapply(myvar, do_kruskal2)
o2 <- bind_rows(o1)
av2$num <- 1:nrow(av2)

uv4 <- merge(av2, o2[,c(5,2)], by.x='Variable', by.y='var', all=T)
uv4 <- uv4[order(uv4$num),]

write.csv2(uv4, file='CRP_baseline_differences.csv')

#=====================================================================
#REVEAL risk scores
r1 <- treatment_response %>% select(STHNumber, TLCOp, PFTDate, DateofFinalDiagnosis)
r1 <- r1 %>% filter(PFTDate !='shared care')
r1 <- r1 %>% filter(PFTDate !='no pfts')
r1 <- r1 %>% filter(PFTDate !='no pft')
r1 <- r1 %>% filter(PFTDate !='no letter')
r1 <- r1 %>% filter(PFTDate !='no  pft')

summary(as.factor(r1$PFTDate))
r1$PFTDate <- as.numeric(r1$PFTDate)

r1$PFTDate <- as.Date(r1$PFTDate, origin = "1899-12-30")
r1$tdiff <- abs(r1$DateofFinalDiagnosis - r1$PFTDate)

r1$STHNumber <- as.factor(r1$STHNumber)
r1$tdiff <- as.numeric(r1$tdiff)
#only select closest recorded value

r2 <- r1 %>% group_by(STHNumber) %>% filter(tdiff == min(tdiff))
#worked

reveal <- mri_correct %>% select(sth_id, final_primary_ph_diagnosis, date_diagnosis, sex, age, pah_subcategory, survival_time, overall_death, pvr, ra_mean, arterial_systolic, heart_rate, who_functional_class, reveal_score, reveal_score_lite, pericardial_effusion, date_iswt, walking_distance, bnp_date, bnp, egfr, who_functional_class)
reveal <- merge(reveal, r2, by.x='sth_id', by.y='STHNumber', all=T)
#dataframe is ready to use with CRP and BMI! 

crp_bmi <- mdf %>% select(ID, level, BMI)
rev <- merge(reveal, crp_bmi, by.x='sth_id', by.y='ID')

#split df in 2
r1 <- rev %>% select(sth_id, date_diagnosis, date_iswt, walking_distance)
r2 <- rev %>% select(sth_id, date_diagnosis, bnp_date, bnp)

result1 <- r1 %>%
  mutate(Difference = abs(date_diagnosis - date_iswt)) %>%  # Compute absolute difference
  group_by(sth_id) %>%
  filter(
    (is.na(Date) & n() == 1) |  # Retain if Date is NA but the only row for that ID
      Difference == min(Difference, na.rm = TRUE)  # Retain closest date otherwise
  ) %>%
  select(-Difference)

#repeat for BNPs
result2 <- r2 %>%
  mutate(Difference = abs(date_diagnosis - bnp_date)) %>%  # Compute absolute difference
  group_by(sth_id) %>%
  filter(
    (is.na(Date) & n() == 1) |  # Retain if Date is NA but the only row for that ID
      Difference == min(Difference, na.rm = TRUE)  # Retain closest date otherwise
  ) %>%
  select(-Difference)

#merge the results together
res <- merge(result1, result2, by='sth_id', all=T)

rev <- rev[,-c(3, 17:20)]
rev <- merge(rev, res, by='sth_id', all=T)

#take only unique vals
rev <- unique(rev)
#filter out all non group 1
rev <- rev %>% filter(final_primary_ph_diagnosis == 'PAH')

#find the duplicate IDs
duplicated_rows <- rev %>% filter(duplicated(sth_id) | duplicated(sth_id, fromLast = TRUE))
#when duplicated, select youngest age as then diagnosed
duplicated_rows <- duplicated_rows %>% group_by(sth_id) %>% slice(which.min(age))
#now check if still duplicated
dupl2 <- duplicated_rows %>% filter(duplicated(sth_id) | duplicated(sth_id, fromLast = TRUE))
#nothing is double here! 
#great, can merge again now

rev <- rev %>% filter(!sth_id %in% duplicated_rows$sth_id)
rev <- rbind(rev, duplicated_rows)

rev <- rev[,-c(18:20,23,24, 26,27)]

#now add CRP and BMI to total score

#first start with adding and removing vars to score
rev$reveal_plus_CRP <- NA
rev$hl_crp <-  ifelse(rev$level > 5, 'high', 'low')
rev$reveal_plus_CRP <-  ifelse(rev$hl_crp == 'high', 2, -2)

rev$reveal_plus_BMI <- 0
rev$reveal_plus_BMI <- ifelse(rev$BMI <18.5, 2, rev$reveal_plus_BMI)
rev$reveal_plus_BMI <- ifelse(rev$BMI >=25 & rev$BMI <30, -1, rev$reveal_plus_BMI)
rev$reveal_plus_BMI <- ifelse(rev$BMI >=30, -2, rev$reveal_plus_BMI)

#now get extra cols
rev$rev_bmi <- rev$reveal_plus_BMI + rev$reveal_score
rev$rev_crp <- rev$reveal_plus_CRP + rev$reveal_score
rev$rev_bmicrp <- rev$reveal_plus_BMI + rev$reveal_score + rev$reveal_plus_CRP

rev$sex_score <- ifelse(rev$sex == 'M' & rev$age >60, 2,0)
rev$renal_score <- ifelse(rev$egfr <60, 1, 0)
rev$sbp_score <- ifelse(rev$arterial_systolic <110, 1, 0)
rev$walk_score <- ifelse(rev$walking_distance >=440, -2, -1)
rev$walk_score <-ifelse(rev$walking_distance <320 & rev$walking_distance >=165, 0, rev$walk_score)
rev$walk_score <- ifelse(rev$walking_distance <165, 1, rev$walk_score)
rev$dlco_score <- ifelse(rev$TLCOp <40, 1, 0)
rev$pe_score <- ifelse(rev$pericardial_effusion == 'Yes', 1, 0)
rev$bnp_score <- ifelse(rev$bnp < 50, -2, 0)
rev$bnp_score <- ifelse(rev$bnp < 800 & rev$bnp >=200, 1, rev$bnp_score)
rev$bnp_score <- ifelse(rev$bnp >= 800, 2, rev$bnp_score)


#make cols 0 if nothing recorded
rev <- rev %>% mutate(across(c(sex_score, renal_score, walk_score, sbp_score, dlco_score, pe_score, bnp_score, reveal_plus_CRP), ~replace_na(., 0)))

#also get non-invasive reveal score ready
rev$noninvasive_reveal <- rev$sex_score + rev$renal_score + rev$sbp_score + rev$walk_score + rev$dlco_score + rev$pe_score + rev$bnp_score
rev$noninvasive_reveal_crp <- rev$sex_score + rev$renal_score + rev$sbp_score + rev$walk_score + rev$dlco_score + rev$pe_score + rev$bnp_score + rev$reveal_plus_CRP
rev$noninvasive_reveal_bmi <- rev$sex_score + rev$renal_score + rev$sbp_score + rev$walk_score + rev$dlco_score + rev$pe_score + rev$bnp_score + rev$reveal_plus_CRP + rev$reveal_plus_BMI


#now get AUC plotted
library(survivalROC)
library(timeROC)

#get survival in years
rev$survival_time <- rev$survival_time/365.2422


#now for 1 year survival
roc1 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$reveal_score, predict.time = 1, method = "KM")
roc2 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$rev_crp, predict.time = 1, method = "KM")
roc3 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$rev_bmi, predict.time = 1, method = "KM")
roc4 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$rev_bmicrp, predict.time = 1, method = "KM")
roc5 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$noninvasive_reveal, predict.time = 1, method = "KM")
roc6 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$noninvasive_reveal_crp, predict.time = 1, method = "KM")
roc7 <- survivalROC(Stime = rev$survival_time, status = rev$overall_death, marker = rev$noninvasive_reveal_bmi, predict.time = 1, method = "KM")


#now plot
plot(roc1$FP, roc1$TP, type = "l", col = "yellow", xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC Curves Comparison of REVEAL score 1 year post diagnosis ASPIRE- Group 1 PAH")
lines(roc2$FP, roc2$TP, col = "black")
lines(roc3$FP, roc3$TP, col = "red")
lines(roc4$FP, roc4$TP, col = "green")
lines(roc5$FP, roc5$TP, col = "blue")
lines(roc6$FP, roc6$TP, col = "cyan")
lines(roc7$FP, roc7$TP, col = "magenta")

legend("bottomright", legend = c(paste("REVEAL 2.0, AUC:", round(roc1$AUC, 2)),
                                 paste("REVEAL + CRP, AUC:", round(roc2$AUC, 2)),
                                 paste("REVEAL + BMI, AUC:", round(roc3$AUC, 2)),
                                 paste("REVEAL + BMI & CRP, AUC:", round(roc4$AUC, 2)),
                                 paste("Noninvasive REVEAL, AUC:", round(roc5$AUC, 2)),
                                 paste("Noninvase REVEAL + CRP, AUC:", round(roc6$AUC, 2)),
                                 paste("Noninvasive REVEAL + CRP & BMI, AUC:", round(roc7$AUC, 2))),
       col = c("yellow", "black", "red", "green", "blue", "cyan", "magenta"), lty = 1, cex=0.8)

