#CRP survival analysis censoring at 10 years post diagnosis
#CRP modelling without clustering
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
getwd()

#load packages
library("lubridate")
library("ggsci")
library("magrittr")
library("tidyverse")
library(glmnet)
library(dplyr)
library(survival)
library(survminer)

#import the data
load("C:/filepath.RData")
crp_df <- data.clean.cohort %>% dplyr::select(id, visit, id_cohort, cbt_inflammation_crp_mgpl, cbt_inflammation_scrp_mgpl)
crp_df <- as.data.frame(crp_df)
#===================================================================================
#Mortality script
#===================================================================================
#import the data
centre <- data.clean.cohort %>% dplyr::select('id', 'centre')
centre <- centre %>% filter(centre %in% c('Glasgow', 'Sheffield', 'Great Ormond Street', 'Lincoln', 'Papworth', 'Royal Brompton', 'Royal United Hospital Bath', 'Imperial and Hammersmith', 'Newcastle Freeman', 'Royal Free'))

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/filepath.rds")

mort_df <-  v4_clean_clinical_data_first_visit_15Dec23 %>% dplyr::select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'bs_bmi', 'DOB', 'sub_cause', 'sub_date')
mort_df <- merge(mort_df, centre)
mort_df$DOB <- as.Date(mort_df$DOB)
mort_df$sub_date <- as.Date(mort_df$sub_date)
mort_df$diagnosis_date <- mort_df$DOB + (mort_df$age_diagnosis*365.2422)
mort_df$diagnosis_date <- as.Date(mort_df$diagnosis_date)
diag_date <- mort_df %>% dplyr::select(id, diagnosis_date)
mort_df <- unique(mort_df)
#make sure sub_date is complete and format the data
library(data.table)
mort_df$sub_date <- as.Date(mort_df$sub_date)
mort_df$sub_date2 <- fifelse(is.na(mort_df$sub_date), ymd("2022-07-01"), mort_df$sub_date)
mort_df$surv_time <- (mort_df$sub_date2 - mort_df$diagnosis_date)
mort_df[c('surv_time', 'day')] <- str_split_fixed(mort_df$surv_time, ' ', 2)
mort_df$surv_time <- as.numeric(mort_df$surv_time)
mort_df$surv_time <- mort_df$surv_time/365.2422
mort_df$event <- ifelse(mort_df$sub_cause == 'death', 1, 0)
mort_df$event <- ifelse(is.na(mort_df$event), 0, mort_df$event)
mort_df$diagnosis_verified <- as.factor(mort_df$diagnosis_verified)
mort_df$sex <- as.factor(mort_df$sex)
mort_df <- mort_df %>% filter(!is.na(surv_time))
#now censor the survival times
mort_df$event <- ifelse(mort_df$surv_time >=10, 0, mort_df$event)
mort_df$surv_time <- ifelse(mort_df$surv_time >=10, 10, mort_df$surv_time)
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)
#=========================================================================================
#now prepare the data for further analyses 
#there is also hs-CRP data --> if CRP is NA use that 
crp_df$crp_diff <- crp_df$cbt_inflammation_crp_mgpl - crp_df$cbt_inflammation_scrp_mgpl
summary(crp_df$crp_diff)
#Differences are minimal
#I am aware of bias risk but roughly corresponds
crp_df$cbt_inflammation_crp_mgpl <- ifelse(is.na(crp_df$cbt_inflammation_crp_mgpl), crp_df$cbt_inflammation_scrp_mgpl, crp_df$cbt_inflammation_crp_mgpl)
crp_df <- crp_df %>% filter(!is.na(cbt_inflammation_crp_mgpl))


#also add additional labels
crp_df$above_5 <- ifelse(crp_df$cbt_inflammation_crp_mgpl >5, 'yes', 'no')

#first CRP_high_low
hl_crp <- crp_df %>% 
  group_by(id) %>% 
  slice_min(order_by = visit)
a<- unique(hl_crp$id)
#this extracts the first visit for everyone! 

#now repeat for high/low CRP
summary(hl_crp$visit)

m_hl <- hl_crp %>% filter(visit <2)

#run mortality analysis
mdf <- merge(mort_df, m_hl)
mdf$above_5 <- as.factor(mdf$above_5)

#remove the kids
mdf <- mdf %>% filter(age_diagnosis>=18)

sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)

#also plot with confidence intervals
ggsurvplot(sfit, data=mdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences based on first recorded CRP in PAH', xlim=c(0,10), legend.title='CRP >=5', break.x.by=2.5, legend.labs=c('First recorded CRP <=5', 'First recorded CRP >5'), palette = c('darkblue', 'darkred'))

#now run with just IPAH
mdf <- mdf %>% filter(diagnosis_verified == 'IPAH')
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
ggsurvplot(sfit, data=mdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences based on first recorded CRP in IPAH', xlim=c(0,10), legend.title='CRP >5', break.x.by=2.5, legend.labs=c('First recorded CRP <=5', 'First recorded CRP >5'), palette = c('darkred','darkblue'))





#now repeat CRP analysis with exclusion of initial high CRPs
excl_crp <- crp_df
#get median CRP
median(excl_crp$cbt_inflammation_crp_mgpl)
IQR(excl_crp$cbt_inflammation_crp_mgpl)
#be stringent, set IQR of 4 as threshold
iqr_4 <- 4*IQR(excl_crp$cbt_inflammation_crp_mgpl) + median(excl_crp$cbt_inflammation_crp_mgpl)

excl_crp <- excl_crp %>% filter(cbt_inflammation_crp_mgpl <iqr_4)

excl_crp <- excl_crp %>% 
  group_by(id) %>% 
  slice_min(order_by = visit)
a<- unique(hl_crp$id)
#this extracts the first visit for everyone! 

#do not take this further than 3rd visit
excl_crp <- excl_crp %>% filter(visit <4)
#1054 patients still included

mdf <- merge(mort_df, excl_crp)
mdf$above_5 <- as.factor(mdf$above_5)
mdf <- mdf %>% filter(age_diagnosis>=18)
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + bs_bmi + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)


#also plot without confidence intervals
ggsurvplot(sfit, data=mdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences based on CRP, after removal of outliers (*4 IQR from median)', xlim=c(0,10), legend.title='CRP >=5', break.x.by=2.5, legend.labs=c('First recorded CRP <=5', 'First recorded CRP >5'), palette = c('darkblue', 'darkred'))

#======================================
#EdB continued here on 17-01-2024
#compare comorbidity between groups
cleaned_comorbidity_dataV2_20Dec2023 <- readRDS("C:/filepath.rds")

com <- merge(m_hl[,c(1,7)], cleaned_comorbidity_dataV2_20Dec2023)
com$above_5 <- as.factor(com$above_5)
com$cfh_comorbid_disease_diagnosis <- as.factor(com$cfh_comorbid_disease_diagnosis)
com$disease_class <- as.factor(com$disease_class)
summary(as.factor(m_hl$above_5))

#get percentage of disease per group
com$n <- 1
com2 <- com %>% group_by(above_5, cfh_comorbid_disease_diagnosis) %>% dplyr::summarise(number = n())
com2$percentage <- ifelse(com2$above_5 == 'yes', (com2$number/381)*100, (com2$number/668)*100)

com3 <- com2[,c(1,2,4)] %>% pivot_wider(names_from = above_5, values_from = percentage)

#get Charlson Comorbidity index
#based on: https://www.mdcalc.com/calc/3917/charlson-comorbidity-index-cci
mi_score <- c('myocardial infarction', 'PCI', 'CABG')
chf_score <- c('cv_heart_failure')
#excluding pulmonary/bronchial aneurysms here
pvd_score <- c('peripheral_vascular_disease', 'abdominal_aortic_aneurysm', 'type_A_aneurysm', 'aortic aneurysm')
cva_score <- c('cva', 'tia')
#no one with dementia in cohort
copd_score <- c('copd')
peptic_ulcer_score <- c('duodenal ulcer', 'gastric ulcer', 'oesophageal ulcer')
liver_score <- c('hepatological')
ctd_score <- c('connective tissue disease', 'scleroderma', 'undifferentiated connective tissue disease', 'antiphospholipid_syndrome', 'overlap syndrome', 'crest syndrome')
hemiplegia_score <- c('hemiparesis')
#CKD score if on dialysis or creatinine >265.26
#for now just give everyone with 'CKD' diagnosis this score
ckd_score <- c('CKD')
solid_tumor_score <- c('solid_malignancy')
leukemia__lymphoma_score <- c('haematological_malignant_and_premalignant')
dm_score <- c('diabetes')
#no AIDS recorded in comorbidities

score_df <- com[,c(1,3,5)]
score_df$score <- 0
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% mi_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$disease_class %in% chf_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% pvd_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% cva_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% copd_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% ctd_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% peptic_ulcer_score, score_df$score+1, score_df$score+0)
#assuming all is mild (lacks data granularity to say more)
score_df$score <- ifelse(score_df$disease_class %in% liver_score, score_df$score+1, score_df$score+0)
#also assuming DM is mild (as end organ damage is difficult to assess easily)
score_df$score <- ifelse(score_df$disease_class %in% dm_score, score_df$score+1, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% hemiplegia_score, score_df$score+2, score_df$score+0)
score_df$score <- ifelse(score_df$cfh_comorbid_disease_diagnosis %in% ckd_score, score_df$score+2, score_df$score+0)
score_df$score <- ifelse(score_df$disease_class %in% solid_tumor_score, score_df$score+2, score_df$score+0)
score_df$score <- ifelse(score_df$disease_class %in% leukemia__lymphoma_score, score_df$score+2, score_df$score+0)
#make sure diseases in one class are not double registrations (e.g. if PCI and MI were both diagnosed)
score_df <- unique(score_df[,c(1,3,4)])
score_df2 <- aggregate(score_df$score, by=list(id=score_df$id), FUN=sum)

clin123 <- v4_clean_clinical_data_first_visit_15Dec23
clin4 <- merge(score_df2, clin123, by='id', all=T)
clin4$x <- ifelse(is.na(clin4$x), 0, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis <50, clin4$x + 0, clin4$x +1)
clin4$x <- ifelse(clin4$age_diagnosis >=60 & clin4$age_diagnosis <= 69, clin4$x + 1, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis >=70 & clin4$age_diagnosis <= 79, clin4$x + 2, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis >=80, clin4$x + 3, clin4$x)

clin4$CCI_adj_score <- clin4$x
clincci <- clin4 %>% select(id, CCI_adj_score)
#=============================
#add cox-PH with morbidity score
mdf <- merge(mort_df, m_hl)
mdf <- mdf %>% filter(age_diagnosis >=18)
mdf <- merge(mdf, clincci, by='id')
names(mdf)[[21]] <- 'CCI_score'
mdf2 <- mdf %>% filter(diagnosis_verified== 'IPAH')

sobj <- Surv(mdf$surv_time, mdf$event, type='right')

#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + bs_bmi + CCI_score, data=mdf)
ggforest(cox, data=mdf)

sobj <- Surv(mdf2$surv_time, mdf2$event, type='right')

cox <- coxph(sobj ~ above_5 + sex + bs_bmi + CCI_score, data=mdf2)
ggforest(cox, data=mdf2)

#use obesity classes
mdf$weight <- ifelse(mdf$bs_bmi <18.5, 'underweight', NA)
mdf$weight <- ifelse(mdf$bs_bmi >=18.5 & mdf$bs_bmi <25, 'normal weight', mdf$weight)
mdf$weight <- ifelse(mdf$bs_bmi >=25 & mdf$bs_bmi <30, 'pre-obesity', mdf$weight)
mdf$weight <- ifelse(mdf$bs_bmi >=30 &  mdf$bs_bmi <35, 'obesity class-I', mdf$weight)
mdf$weight <- ifelse(mdf$bs_bmi >=35 &  mdf$bs_bmi <40, 'obesity class-II', mdf$weight)
mdf$weight <- ifelse(mdf$bs_bmi >=40, 'obesity class-III', mdf$weight)
mdf$weight <- as.factor(mdf$weight)

sobj <- Surv(mdf$surv_time, mdf$event, type='right')

#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + weight + CCI_score, data=mdf)
ggforest(cox, data=mdf)
