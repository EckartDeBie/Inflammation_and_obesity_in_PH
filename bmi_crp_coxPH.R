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
definitive.clusterings.3k.PAM.labels.3March2021 <- read.csv("~/Not PhD/My publications/Data new analyses March 2021/definitive clusterings 3k PAM labels 3March2021.csv", sep=";")
load("C:/Users/Gebruiker/Downloads/data_clean (3).RData")
crp_df <- data.clean.cohort %>% select(id, visit, id_cohort, cbt_inflammation_crp_mgpl, cbt_inflammation_scrp_mgpl)
crp_df <- as.data.frame(crp_df)

#===================================================================================
#Mortality script
#===================================================================================
#import the data
centre <- data.clean.cohort %>% select('id', 'centre')
centre <- centre %>% filter(centre %in% c('Glasgow', 'Sheffield', 'Great Ormond Street', 'Lincoln', 'Papworth', 'Royal Brompton', 'Royal United Hospital Bath', 'Imperial and Hammersmith', 'Newcastle Freeman', 'Royal Free'))

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")

mort_df <-  v4_clean_clinical_data_first_visit_15Dec23 %>% select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'bs_bmi', 'DOB', 'sub_cause', 'sub_date')
mort_df <- merge(mort_df, centre)
mort_df$DOB <- as.Date(mort_df$DOB)
mort_df$sub_date <- as.Date(mort_df$sub_date)
mort_df$diagnosis_date <- mort_df$DOB + (mort_df$age_diagnosis*365.2422)
mort_df$diagnosis_date <- as.Date(mort_df$diagnosis_date)
diag_date <- mort_df %>% select(id, diagnosis_date)
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
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)
#=========================================================================================
#impute with HS-crp
crp_df$cbt_inflammation_crp_mgpl <- ifelse(is.na(crp_df$cbt_inflammation_crp_mgpl), crp_df$cbt_inflammation_scrp_mgpl, crp_df$cbt_inflammation_crp_mgpl)
crp_df <- crp_df %>% filter(!is.na(cbt_inflammation_crp_mgpl))


#also add additional labels
crp_df$above_5 <- ifelse(crp_df$cbt_inflammation_crp_mgpl >5, 'yes', 'no')


#first CRP_high_low
hl_crp <- crp_df %>% 
  group_by(id) %>% 
  slice_min(order_by = visit)

hl_lab <- hl_crp %>% filter(visit <2) %>% select(id, above_5)

#make BMI factorial
hl_mort <- merge(mort_df, hl_lab)
bmi_df <- hl_mort
bmi_df$weight <- ifelse(bmi_df$bs_bmi <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=18.5 & bmi_df$bs_bmi <25, 'normal weight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=25 & bmi_df$bs_bmi <30, 'pre-obesity', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=30 &  bmi_df$bs_bmi <35, 'obesity class-I', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=35 &  bmi_df$bs_bmi <40, 'obesity class-II', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=40, 'obesity class-III', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)
df_save <- bmi_df %>% select(id, weight)
saveRDS(df_save,'BMI_category_per_patient.rds')

hl_mort <- bmi_df
hl_mort$above_5 <- as.factor(hl_mort$above_5)

mdf <- hl_mort
mdf <- mdf %>% filter(!is.na(weight))

CCI_score_per_patient <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CCI_score_per_patient.rds")

mdf <- merge(mdf, CCI_score_per_patient, by='id')

sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + weight + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)

#include just BMI
cox <- coxph(sobj ~  sex + age_diagnosis + weight, data=mdf)
ggforest(cox, data=mdf)

#repeat with BMI as numeric variable
cox <- coxph(sobj ~  sex + age_diagnosis + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

#limit survival analysis to 5 years
mdf2 <- mdf
mdf2$event <- ifelse(mdf2$surv_time >=5, 0, mdf2$event)
mdf2$surv_time <- ifelse(mdf2$surv_time >=5, 5, mdf2$surv_time)
levels(mdf2$weight)
#relevel
mdf2$weight2 <- mdf2$weight
mdf2$weight2 <- factor(mdf2$weight2, levels = c("underweight", "normal weight", "pre-obesity", "obesity class-I", "obesity class-II", "obesity class-III"))
levels(mdf2$weight2)

#repeat cox-PH
cox <- coxph(sobj ~  sex + age_diagnosis + weight, data=mdf2)
ggforest(cox, data=mdf2)

#scale and repeatCoxPH

hd <- v4_clean_clinical_data_first_visit_15Dec23 %>% select(id, hb_pvr_calc, hb_cardiac_index_value_1, hb_rap_m, hb_cardiac_output_value_1, hb_pawp_m)
mdf123 <- merge(mdf2, hd, by='id', all=T)

mdf123$pvr <- mdf123$hb_pvr_calc/80

mdf123$scale_age <- scale(mdf123$age_diagnosis)
mdf123$scale_cci <- scale(mdf123$CCI_adj_score)
mdf123$scale_pvr <- scale(mdf123$pvr)
mdf123$scale_rap <- scale(mdf123$hb_rap_m)
#filter out kids
mdf123 <- mdf123 %>% filter(age_diagnosis >=18)

#n=792 remain
sobj <- Surv(mdf123$surv_time, mdf123$event, type='right')

#now do coxPH #and add in smoking Hx
query_smoking_dump_10_7_23 <- read.delim("~/PhD/Projects/CRP ~ survival and BMI/query_smoking_dump_10_7_23.txt")
smoke <- query_smoking_dump_10_7_23
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_history == 'NULL', 'no', smoke$cfh_smoking_history)
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_history == 'UNK', NA, smoke$cfh_smoking_history)
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_current == 'yes', 'Current smoker', smoke$cfh_smoking_history)
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_history == 'yes', 'Smoking history', smoke$cfh_smoking_history)

smoke$`Smoking history` <- as.factor(smoke$cfh_smoking_history)
#merge with CoxPH data
mdf1234 <- merge(mdf123, smoke[,c(2,6)], by.x='id', by.y='study_subject_oid')
#losing 12 people but that's fine I think....

sobj <- Surv(mdf1234$surv_time, mdf1234$event, type='right')

levels(mdf1234$`Smoking history`)
mdf1234$`Smoking history` <- relevel(mdf1234$`Smoking history`, ref=3)

cox <- coxph(sobj ~  weight + sex + `Smoking history`+ scale_age + scale_cci + scale_pvr + scale_rap, data=mdf1234)
ggforest(cox, data=mdf1234)
sc <- summary(cox)
sc123 <- as.data.frame(sc$coefficients)
saveRDS(sc123, file='SCALED_coef_IPAH_cohort_BMI_groups.rds')

library(contsurvplot)

#=================================================================================
#get a clinical differences table (like in ASPIRE)
#be consistent and take all pts
v4 <- v4_clean_clinical_data_first_visit_15Dec23

v4$weight <- ifelse(v4$bs_bmi <18.5, 'underweight', NA)
v4$weight <- ifelse(v4$bs_bmi >=18.5 & v4$bs_bmi <25, 'normal weight', v4$weight)
v4$weight <- ifelse(v4$bs_bmi >=25 & v4$bs_bmi <30, 'pre-obesity', v4$weight)
v4$weight <- ifelse(v4$bs_bmi >=30 &  v4$bs_bmi <35, 'obesity class-I', v4$weight)
v4$weight <- ifelse(v4$bs_bmi >=35 &  v4$bs_bmi <40, 'obesity class-II', v4$weight)
v4$weight <- ifelse(v4$bs_bmi >=40, 'obesity class-III', v4$weight)
v4$weight <- as.factor(v4$weight)

library(Publish)
uv1 <- univariateTable(weight ~ Gender + AgeAtDiag + FinalPrimaryPHDiagnosis + BMI + level + PVR + pvri + tpr + CardiacIndex + CardiacOutput + mRAP + mPAP + PAWP + FVCp + eGFR + BaseISWD, data=mdf, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE, summary.format = 'median(x) [iqr(x)]')
uv2 <- summary(uv1)

uv3 <- uv2[,c(1,2,8,3,7,4:6, 9, 10)]

#========================================================================================================================
#get baseline data!
#=================================================================================
#now repeat baseline table!
library(Publish)
library(broom)
b <- v4

b <- b %>% filter(!is.na(weight))
b <- b %>% filter(age_diagnosis >=18)

b$tpr <- b$hb_pap_m/b$hb_cardiac_output_value_1
b$pvri <- (b$hb_pap_m - b$hb_pawp_m)/b$hb_cardiac_index_value_1

a1 <- univariateTable(weight ~ diagnosis_verified + sex + age_diagnosis + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + tpr + pvri + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc + functional_class + cbt_card_ntprobnp_ngpl, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b)
a2 <- summary(a1)


#write.csv2(a2, file='univar_table_log_crp_clusters.csv')

do_chisq <- function(data, myvar, group_var) {
  tbl <- data[[myvar]]
  grp <- data[[group_var]]
  chi_output <- chisq.test(table(tbl, grp))
  chi_output <- tidy(chi_output)
  chi_output$var <- myvar
  return(chi_output)
}

myvar_ch <- c('diagnosis_verified', 'sex', 'hv_vasodilator_responder', 'functional_class')


group_var = "weight"  
output <- lapply(myvar_ch, function(x) do_chisq(b, x, group_var))
combined_output1 <- bind_rows(output)


chisq.test(b$diagnosis_verified, b$weight)
#0.005064
chisq.test(b$hv_vasodilator_responder, b$weight)
#0.5933

combined_output1$p.value <- ifelse(combined_output1$var == 'diagnosis_verified', 0.005064, combined_output1$p.value)
combined_output1$p.value <- ifelse(combined_output1$var == 'hv_vasodilator_responder', 0.5933, combined_output1$p.value)


#make function for kruskal test

do_mw <- function(myvar) {
  mw_output <- kruskal.test(b[[myvar]] ~ b$weight)
  mw_output <- (tidy(mw_output))
  mw_output$var <- myvar
  return(mw_output)
}

do_mw('age_diagnosis')

#now do for all
myvar_mw <- c('age_diagnosis', 'cbt_thyr_tsh_mupl', 'tpr', 'pvri', 'cbt_thyr_freet4_pmolpl', 'egfr_mdrd', 'cbt_haem_platelets_x10e9pl', 'cbt_haem_hb_gpl', 'cfe_rest_spo2', 'cfe_heart_rate', 'hb_pawp_m', 'hb_pap_d',  'hb_pap_m', 'hb_pvr_calc', 'hb_rap_m', 'hb_cardiac_output_value_1', 'hb_cardiac_index_value_1', 'lf_fev1_pc', 'lf_fvc_pc', 'lf_kco_pc', 'cbt_card_ntprobnp_ngpl')
myvar = myvar_mw
output <- lapply(myvar, do_mw)
combined_output2 <- bind_rows(output)

pvals <- rbind(combined_output1, combined_output2)

#quick repeat with just corridor walk
b2 <- b %>% filter(ep_1_type_6mwt == 'corridor')
c1 <- univariateTable(weight ~ ep_1_distance_meters, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b2)
c2 <- summary(c1)

kruskal.test(b2$ep_1_distance_meters ~ b2$weight)
#pval =  2.2e-16

pvals[26,] <- NA
pvals$p.value <- ifelse(is.na(pvals$var),  2.2e-16, pvals$p.value)
pvals$var <- ifelse(is.na(pvals$var), 'ep_1_distance_meters', pvals$var)

#get pval for CCI scores (from plot below)
pvals[27,] <- NA
pvals$p.value <- ifelse(is.na(pvals$var), 1e-11, pvals$p.value)
pvals$var <- ifelse(is.na(pvals$var), 'CCI_score', pvals$var)

pvals$fdr <- p.adjust(pvals$p.value, method='fdr')
pvals$fdr_round <- round(pvals$fdr, digits=10)


#get a CCI table
CCI_score_per_patient_complete <- readRDS("C:/Users/location/CCI_score_per_patient_complete.rds")
b3 <- merge(b, CCI_score_per_patient_complete, by='id')

d1 <- univariateTable(weight ~ CCI_adj_score, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b3)
d2 <- summary(d1)

ac <- rbind(a2, d2)


p8 <- ggplot(b3, aes(x=weight, y=CCI_adj_score, fill=weight)) + 
geom_boxplot(outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in Charlson Comorbidity Index between CRP groups') +
  ylab('CCI') +
  xlab('CRP group') +
  theme_bw() +
  stat_compare_means(method = 'kruskal.test') +
  scale_x_discrete(labels= labs)

#now merge these data with the univariate table

pv2 <- pvals %>% select(var, p.value, fdr_round)

ac$increasing_numbers <- seq(1, nrow(ac))


a3 <- merge(ac, pv2, by.x='Variable', by.y='var', all=T)

a3 <- a3[order(a3$increasing_numbers), ]

write.csv2(a3, file='univar_table_baseline_weight_groups.csv')

#===================================================================================


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

#============================================================================


#repeat with BMI as numeric variable
cox <- coxph(sobj ~  sex + age_diagnosis + bs_bmi, data=mdf2)
ggforest(cox, data=mdf2)

s1 <- Surv(mdf2$surv_time, mdf2$event, type='right')
s2 <- survfit(s1 ~ weight, data=mdf2)

ggsurvplot(s2, data=mdf2, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = FALSE, title='survival differences based on BMI', xlim=c(0,5), legend.title='Group', break.x.by=0.5, legend.labs=c('underweight', 'normal weight', 'pre-obesity', 'obesity class-I', 'obesity class-II', 'obesity class-III'))


#now remove kids and only take IPAH
mdf2 <- mdf2 %>% filter(age_diagnosis > 18 & diagnosis_verified == 'IPAH')

s1 <- Surv(mdf2$surv_time, mdf2$event, type='right')
s2 <- survfit(s1 ~ weight, data=mdf2)

ggsurvplot(s2, data=mdf2, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = FALSE, title='survival differences based on BMI', xlim=c(0,5), legend.title='Group', break.x.by=0.5, legend.labs=c('underweight', 'normal weight', 'pre-obesity', 'obesity class-I', 'obesity class-II', 'obesity class-III'))

#repeat cox-PH
cox <- coxph(s1 ~  sex + age_diagnosis + weight, data=mdf2)
ggforest(cox, data=mdf2)

#repeat with BMI as numeric variable
cox <- coxph(s1 ~  sex + age_diagnosis + bs_bmi, data=mdf2)
ggforest(cox, data=mdf2)

cox <- coxph(s1 ~  sex + age_diagnosis + bs_bmi + CCI_adj_score, data=mdf2)
ggforest(cox, data=mdf2)

#include haemodynamics
hd <- v4_clean_clinical_data_first_visit_15Dec23 %>% select(id, hb_pvr_calc, hb_cardiac_index_value_1, hb_rap_m, hb_cardiac_output_value_1, hb_pawp_m)
mdf4 <- merge(mdf2, hd, by='id')
mdf4$pvri <- mdf4$hb_pvr_calc/mdf4$hb_cardiac_index_value_1

s1 <- Surv(mdf4$surv_time, mdf4$event, type='right')

cox <- coxph(s1 ~  sex + age_diagnosis + bs_bmi + CCI_adj_score + pvri, data=mdf4)
ggforest(cox, data=mdf4)


#repeat the analysis with just IPAH and >18 at diagnosis
mdf <- mdf %>% filter(age_diagnosis > 18 & diagnosis_verified == 'IPAH')
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + weight + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)

#compare with clusters
md2 <- merge(hl_mort, data.clean.cohort, by='id')
md2 <- md2 %>% select(id, id_cohort, above_5)
md2 <- unique(md2)

labs_clust <- merge(md2, definitive.clusterings.3k.PAM.labels.3March2021, by.x='id_cohort', by.y='X')
library(Publish)
labs_clust$definitive_labels <- as.factor(labs_clust$definitive_labels)
a1 <- univariateTable(definitive_labels ~ above_5, data=labs_clust, column.percent = TRUE, show.totals = TRUE, compare.groups = TRUE)
a2 <- summary(a1)

cox <- coxph(sobj ~ weight + age_diagnosis + sex, data=mdf)
ggforest(cox, data=mdf)

#assess treatment effects
treatment_effects_file_v1_March2024 <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/treatment_effects_file_v1_March2024.rds")

wl <- mdf %>% select(id, weight)
wl <- merge(wl, treatment_effects_file_v1_March2024, by='id')

summary(aov(wl$walk_diff_meters ~ wl$weight))
summary(aov(wl$fc_diff ~ wl$weight))



