#linear modelling of variables based on CRP/BMI
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

#import the data
load("C:/Users/Gebruiker/Downloads/data_clean (3).RData")
df_c <- data.clean.cohort %>% select(id, visit, id_cohort, cbt_inflammation_crp_mgpl, cbt_inflammation_scrp_mgpl)
df_c <- as.data.frame(df_c)

#there is also hsCRP data --> if CRP is NA use that 
df_c$crp_diff <- df_c$cbt_inflammation_crp_mgpl - df_c$cbt_inflammation_scrp_mgpl
summary(df_c$crp_diff)
#Differences are minimal
#I am aware of bias risk but roughly corresponds
df_c$cbt_inflammation_crp_mgpl <- ifelse(is.na(df_c$cbt_inflammation_crp_mgpl), df_c$cbt_inflammation_scrp_mgpl, df_c$cbt_inflammation_crp_mgpl)
df_c <- df_c %>% filter(!is.na(cbt_inflammation_crp_mgpl))

crp_dat <- df_c
crp_for_later <- crp_dat

centre <- data.clean.cohort %>% select('id', 'centre')
centre <- centre %>% filter(centre %in% c('Glasgow', 'Sheffield', 'Great Ormond Street', 'Lincoln', 'Papworth', 'Royal Brompton', 'Royal United Hospital Bath', 'Imperial and Hammersmith', 'Newcastle Freeman', 'Royal Free'))

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")
CCI_score_per_patient_complete <- readRDS("C:/Users/location/CCI_score_per_patient_complete.rds")

#========================================================
#repeat Lms with just adults --> generate the data-frames
clin <- merge(v4_clean_clinical_data_first_visit_15Dec23, CCI_score_per_patient_complete, by='id')
clin <- clin %>% filter(age_diagnosis >=18)


hl_crp <- crp_for_later %>% 
  group_by(id) %>% 
  slice_min(order_by = visit)

#only select CRP from the first 3 visits to maximise results
hl_crp <- hl_crp %>% filter(visit <3)
hl_crp <- hl_crp %>% select(id, cbt_inflammation_crp_mgpl)
names(hl_crp) <- c('id', 'pref_crp')

clin <- merge(hl_crp, clin)

#subset the clinical file
clin1 <- clin
clusters_unscaled_crp_time_pam_2k <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/clusters_unscaled_crp_time_pam_2k.rds")
clin2 <- clin
clin2 <- clin2 %>% filter(!id %in% clusters_unscaled_crp_time_pam_2k$name)

plot(clin1$cbt_inflammation_crp_mgpl, clin$bs_bmi)
cor.test(log(clin1$cbt_inflammation_crp_mgpl), log(clin1$bs_bmi))
#================================================================================
#first plot age and BMI to get interaction

dob <- v4_clean_clinical_data_first_visit_15Dec23 %>% select(id, DOB)

bm_df <- data.clean.cohort %>% select(id, visit, bs_bmi, bs_date)
mis_df <- data.clean.cohort %>% select(id, visit, bs_bmi, bs_date)
summary(mis_df$bs_bmi)
#1240/4629 repeat BMIs missing --> 27%
bm_df <- merge(bm_df, dob, by='id')
bm_df$bs_date <- as.Date(bm_df$bs_date)
bm_df$DOB <- as.Date(bm_df$DOB)
bm_df$bmi_age <- as.numeric(bm_df$bs_date - bm_df$DOB)/365.2422

#clean up the df
bm_df <- bm_df %>% filter(!is.na(bs_bmi))
hist(bm_df$bs_bmi)
#exclude BMI >1000
library(naniar)
bm_df <- bm_df %>% replace_with_na_at(.vars = c('bs_bmi'), condition = ~.x >250)
#now looks roughly normal with some outliers

ggplot(bm_df, aes(bmi_age, bs_bmi)) +
  geom_point() +
  geom_smooth() +
  ggtitle('BMI by age in PAH') +
  labs(y='BMI (kg/m^2)', x='age (years)')

ggplot(bm_df, aes(visit, bs_bmi)) +
  geom_point() +
  geom_smooth() +
  ggtitle('BMI by visit') +
  labs(y='BMI (kg/m^2)', x='Study visit')

#=========================================================================================
#import mortality data
mort_df <-  v4_clean_clinical_data_first_visit_15Dec23 %>% select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'DOB', 'sub_cause', 'sub_date')
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
#=================================================================================
#see if BMI changes for the survivers vs non-survivers
mortality_labs <- mort_df %>% select(id, sub_cause)
#remove everyone who is not dead or NA (LTx Loss to FU etc all removed)
mortality_labs <- mortality_labs %>% filter(is.na(sub_cause) | sub_cause == 'death')
mortality_labs$sub_cause <- ifelse(is.na(mortality_labs$sub_cause), 'alive', 'dead')
mortality_labs$sub_cause <- as.factor(mortality_labs$sub_cause)

bm_df <- merge(bm_df, mortality_labs, by='id')

ggplot(bm_df, aes(x=bmi_age, y=bs_bmi, fill=sub_cause)) +
  geom_point() +
  geom_smooth() +
  ggtitle('BMI by age in PAH') +
  labs(y='BMI (kg/m^2)', x='age (years)')
#repeat with removal of kids

bm_df_adult <- bm_df %>% filter(bmi_age >18)
ggplot(bm_df_adult, aes(x=bmi_age, y=bs_bmi, fill=sub_cause)) +
  geom_point() +
  geom_smooth() +
  ggtitle('BMI by age in PAH') +
  labs(y='BMI (kg/m^2)', x='age (years)')

ggplot(bm_df, aes(x=visit, y=bs_bmi, fill=sub_cause)) +
  geom_point() +
  geom_smooth() +
  ggtitle('BMI by visit') +
  labs(y='BMI (kg/m^2)', x='Study visit')

#=======================================================================
#now fit some LMs
library(ggResidpanel)
library(broom)

#write a formula
do_lm_log = function(myvar) {
  
#1) now model myvar
#only take complete data to model the same df
df <- clin1 %>% dplyr::select(id, myvar, age_diagnosis, sex, bs_bmi, pref_crp, CCI_adj_score)
df <- na.omit(df)
  
fit_0 = lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex))
fit_1 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$bs_bmi))
fit_2 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp))
fit_3 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$bs_bmi) + log(df$pref_crp))
fit_4 <- lm(log(df[[myvar]]) ~ + df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp)*log(df$bs_bmi))
fit_5 <- lm(log(df[[myvar]]) ~ + df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp) + log(df$bs_bmi) + df$CCI_adj_score)


plot(resid_panel(fit_0,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_1,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_2,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_3,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_4,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_5,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))


print(summary(fit_0))
print(summary(fit_1))
print(summary(fit_2))
print(summary(fit_3))
print(summary(fit_4))
print(summary(fit_5))

print(anova(fit_1, fit_0))
print(anova(fit_2, fit_1))
print(anova(fit_3, fit_2))
print(anova(fit_4, fit_3))
print(anova(fit_3, fit_5))

print(BIC(fit_0, fit_1, fit_2, fit_3, fit_4, fit_5))

fit_5_summary <- summary(fit_5)
tidy_summary <- tidy(fit_5_summary)

tidy_summary$variable_tested <- myvar

return(tidy_summary)

}

do_lm_no_log = function(myvar) {
  
  #1) now model myvar
  #only take complete data to model the same df
  df <- clin1 %>% dplyr::select(id, myvar, age_diagnosis, sex, bs_bmi, pref_crp, CCI_adj_score)
  df <- na.omit(df)
  
fit_0 = lm(df[[myvar]]  ~ df$age_diagnosis + as.factor(df$sex))
fit_1 <- lm(df[[myvar]]  ~ df$age_diagnosis + as.factor(df$sex) + log(df$bs_bmi))
fit_2 <- lm(df[[myvar]]  ~ df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp))
fit_3 <- lm(df[[myvar]]  ~ df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp) + log(df$bs_bmi))
fit_4 <- lm(df[[myvar]] ~ df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp)*log(df$bs_bmi))
fit_5 <- lm(df[[myvar]] ~ + df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp) + log(df$bs_bmi) + df$CCI_adj_score)

  
  plot(resid_panel(fit_0,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_1,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_2,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_3,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_4,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_5,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  
  
  print(summary(fit_0))
  print(summary(fit_1))
  print(summary(fit_2))
  print(summary(fit_3))
  print(summary(fit_4))
  print(summary(fit_5))
  
  
  print(anova(fit_1, fit_0))
  print(anova(fit_2, fit_1))
  print(anova(fit_3, fit_2))
  print(anova(fit_4, fit_3))
  print(anova(fit_3, fit_5))
  
 print(BIC(fit_0, fit_1, fit_2, fit_3, fit_4, fit_5))
 
 fit_5_summary <- summary(fit_5)
 tidy_summary <- tidy(fit_5_summary)

 tidy_summary$variable_tested <- myvar
 
 return(tidy_summary)
}

#log-transformation of CO and CI makes it look better

log_vars <- c('hb_pawp_m', 'hb_pvr_calc', 'hb_cardiac_index_value_1', 'hb_cardiac_output_value_1', 'hb_pap_m', 'cbt_card_ntprobnp_ngpl')
nolog_vars <- c('egfr_mdrd', 'hb_rap_m', 'lf_fvc_pc')

#now do LM for Wedge and RAP

#do_lm_log('hb_pawp_m')
#do_lm_no_log('hb_rap_m')
#do_lm_log('hb_pvr_calc')
#do_lm_no_log('egfr_mdrd')
#do_lm_no_log('lf_fvc_pc')

res <- lapply(log_vars, do_lm_log)
res2 <- do.call(rbind, res)

resnl <- lapply(nolog_vars, do_lm_no_log)
resnl2 <- do.call(rbind, resnl)

lm_results <- rbind(res2, resnl2)

clin1a <- clin1
clin1 <- clin1 %>% filter(ep_1_type_6mwt == 'corridor')
abc <- do_lm_no_log('ep_1_distance_meters')
clin1 <- clin1a

lm_results <- rbind(lm_results, abc)

do_lm_crp = function(myvar) {
  
df <- clin1 %>% dplyr::select(id, age_diagnosis, sex, bs_bmi, pref_crp, CCI_adj_score)
df <- na.omit(df)

fit_0 = lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex))
fit_1 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$bs_bmi))
fit_3 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$bs_bmi) + df$CCI_adj_score)


plot(resid_panel(fit_0,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_1,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
plot(resid_panel(fit_3,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))



print(summary(fit_0))
print(summary(fit_1))
print(summary(fit_3))



print(BIC(fit_0, fit_1, fit_3))

fit_3_summary <- summary(fit_3)
tidy_summary <- tidy(fit_3_summary)

tidy_summary$variable_tested <- myvar

return(tidy_summary)

}


do_lm_bmi = function(myvar) {
  
  df <- clin1 %>% dplyr::select(id, age_diagnosis, sex, bs_bmi, pref_crp, CCI_adj_score)
  df <- na.omit(df)
  
  fit_0 = lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex))
  fit_2 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp))
  fit_4 <- lm(log(df[[myvar]])  ~ df$age_diagnosis + as.factor(df$sex) + log(df$pref_crp) + df$CCI_adj_score)
  
  
  plot(resid_panel(fit_0,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_2,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  plot(resid_panel(fit_4,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE))
  
  
  
  print(summary(fit_0))
  print(summary(fit_2))
  print(summary(fit_4))
  
  
  print(BIC(fit_0, fit_2, fit_4))
  
  fit_4_summary <- summary(fit_4)
  tidy_summary <- tidy(fit_4_summary)
  
  tidy_summary$variable_tested <- myvar
  
  return(tidy_summary)
  
}

resbmi <- do_lm_bmi('bs_bmi')
rescrp <- do_lm_crp('pref_crp')

resshort <- rbind(resbmi, rescrp)
lm_results <- rbind(lm_results, resshort)
write_rds(lm_results, 'results_linear_models_cohort.rds')

lm_res_adj <- lm_results

sd_crp <- sd(log(clin1$pref_crp))
sd_bmi <- clin1 %>% filter(!is.na(bs_bmi))
sd_bmi <- sd(log(sd_bmi$bs_bmi))
sd_pawp <- clin1 %>% filter(!is.na(hb_pawp_m))
sd_pawp <- sd(log(sd_pawp$hb_pawp_m))
sd_pvr <- clin1 %>% filter(!is.na(hb_pvr_calc))
sd_pvr <- sd(log(sd_pvr$hb_pvr_calc))
sd_pap <- clin1 %>% filter(!is.na(hb_pap_m))
sd_pap <- sd(log(sd_pap$hb_pap_m))
sd_bnp <- clin1 %>% filter(!is.na(cbt_card_ntprobnp_ngpl))
sd_bnp <- sd(log(sd_bnp$cbt_card_ntprobnp_ngpl))
sd_co <- clin1 %>% filter(!is.na(hb_cardiac_output_value_1))
sd_co <- sd(log(sd_co$hb_cardiac_output_value_1))
sd_ci <- clin1 %>% filter(!is.na(hb_cardiac_index_value_1))
sd_ci <- sd(log(sd_ci$hb_cardiac_index_value_1))

#also get sd data for non-log-transformed data
sd_fvc <- clin1 %>% filter(!is.na(lf_fvc_pc))
sd_fvc <- sd(sd_fvc$lf_fvc_pc)
sd_egfr <- clin1 %>% filter(!is.na(egfr_mdrd))
sd_egfr <- sd(sd_egfr$egfr_mdrd)
sd_rap <- clin1 %>% filter(!is.na(hb_rap_m))
sd_rap <- sd(sd_rap$hb_rap_m)
sd_mwd <- clin1 %>% filter(!is.na(ep_1_distance_meters) & ep_1_type_6mwt == 'corridor')
sd_mwd <- sd(sd_mwd$ep_1_distance_meters)

lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'bs_bmi', lm_res_adj$estimate/sd_bmi, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'pref_crp', lm_res_adj$estimate/sd_crp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'lf_fvc_pc', lm_res_adj$estimate/sd_fvc, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'hb_rap_m', lm_res_adj$estimate/sd_rap, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'egfr_mdrd', lm_res_adj$estimate/sd_egfr, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'hb_pvr_calc', lm_res_adj$estimate/sd_pvr, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'ep_1_distance_meters', lm_res_adj$estimate/sd_mwd, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'hb_pawp_m', lm_res_adj$estimate/sd_pawp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'hb_pap_m', lm_res_adj$estimate/sd_pap, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'cbt_card_ntprobnp_ngpl', lm_res_adj$estimate/sd_bnp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'hb_cardiac_output_value_1', lm_res_adj$estimate/sd_co, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$variable_tested == 'hb_cardiac_index_value_1', lm_res_adj$estimate/sd_ci, lm_res_adj$estimate)


#also adjust standard error
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'bs_bmi', lm_res_adj$std.error/sd_bmi, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'pref_crp', lm_res_adj$std.error/sd_crp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'lf_fvc_pc', lm_res_adj$std.error/sd_fvc, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'hb_rap_m', lm_res_adj$std.error/sd_rap, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'egfr_mdrd', lm_res_adj$std.error/sd_egfr, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'hb_pvr_calc', lm_res_adj$std.error/sd_pvr, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'ep_1_distance_meters', lm_res_adj$std.error/sd_mwd, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'hb_pawp_m', lm_res_adj$std.error/sd_pawp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'hb_pap_m', lm_res_adj$std.error/sd_pap, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'cbt_card_ntprobnp_ngpl', lm_res_adj$std.error/sd_bnp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'hb_cardiac_output_value_1', lm_res_adj$std.error/sd_co, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$variable_tested == 'hb_cardiac_index_value_1', lm_res_adj$std.error/sd_ci, lm_res_adj$std.error)

lm_res_adj$std.error <- lm_res_adj$std.error*1.96


#also adjust for other tested variable in model
lm_res_adj$std.error <- ifelse(lm_res_adj$term == 'log(df$pref_crp)', lm_res_adj$std.error*sd_crp, lm_res_adj$std.error)
lm_res_adj$std.error <- ifelse(lm_res_adj$term == 'log(df$bs_bmi)', lm_res_adj$std.error*sd_bmi, lm_res_adj$std.error)

lm_res_adj$estimate <- ifelse(lm_res_adj$term == 'log(df$pref_crp)', lm_res_adj$estimate*sd_crp, lm_res_adj$estimate)
lm_res_adj$estimate <- ifelse(lm_res_adj$term == 'log(df$bs_bmi)', lm_res_adj$estimate*sd_bmi, lm_res_adj$estimate)


#filter only the variables for plotting
lm_res_adj <- lm_res_adj %>% filter(term %in% c('log(df$pref_crp)', 'log(df$bs_bmi)'))
saveRDS(lm_res_adj, file='prepared_df_for_forest_models_cohort.rds')


#now calculate PVRI
clin1$pvri <- (clin1$hb_pap_m - clin1$hb_pawp_m)/clin1$ci
hist(clin1$pvri)
#clear skew --> log-transform

do_lm_no_log('egfr_mdrd')
do_lm_no_log('lf_fvc_pc')

do_lm_log('pvri')

#also model 6MWD (this needs pre-selecting so do inidivual model)
#clin1 <- clin1 %>% filter(ep_1_type_6mwt == 'corridor')
#do_lm_no_log('ep_1_distance_meters')

#include smoking Hx in 6MWD
query_smoking_dump_10_7_23 <- read.delim("~/PhD/Projects/CRP ~ survival and BMI/query_smoking_dump_10_7_23.txt")
smoke <- query_smoking_dump_10_7_23
smoke2 <- smoke %>% select(study_subject_oid, cfh_smoking_history, cfh_smoking_current)
smoke2 <- unique(smoke2)
smoke2$cfh_smoking_history <- ifelse(smoke2$cfh_smoking_history == 'NULL', 'no', smoke2$cfh_smoking_history)
smoke2$cfh_smoking_history <- ifelse(smoke2$cfh_smoking_history == 'UNK', NA, smoke2$cfh_smoking_history)
smoke2$cfh_smoking_history <- as.factor(smoke2$cfh_smoking_history)

smoke2$cfh_smoking_current <- ifelse(smoke2$cfh_smoking_current == 'NULL', 'no', smoke2$cfh_smoking_current)
smoke2$cfh_smoking_current <- ifelse(smoke2$cfh_smoking_current == 'UNK', NA, smoke2$cfh_smoking_current)
smoke2$cfh_smoking_current <- as.factor(smoke2$cfh_smoking_current)

clin2 <- merge(clin1, smoke2, by.x='id', by.y='study_subject_oid')

#now model a linear model
fit <- lm(clin2[["ep_1_distance_meters"]]  ~ clin2$age_diagnosis + as.factor(clin2$sex) + log(clin2$pref_crp) + log(clin2$bs_bmi) + as.factor(clin2$ep_1_type_6mwt) + clin2$cfh_smoking_current.y + clin2$cfh_smoking_history.y)
resid_panel(fit,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE)
summary(fit)

BIC(fit)
#===================================================================================
#EdB continued here on 09-02-24
#model CRP based on haemodynamics

haem <- clin1 %>% select(id, bs_bmi, pref_crp, hb_pvr_calc, hb_pawp_m, hb_rap_m, age_diagnosis, sex)
haem <- na.omit(haem)





f1 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi))
f2 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + log(haem$hb_pvr_calc))
f3 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + log(haem$hb_pawp_m))
f4 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + haem$hb_rap_m)
f5 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + log(haem$hb_pvr_calc) + log(haem$hb_pawp_m))
f6 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + haem$hb_rap_m + log(haem$hb_pvr_calc))
f7 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + haem$hb_rap_m + log(haem$hb_pawp_m))
f8 <- lm(log(haem$pref_crp) ~ log(haem$bs_bmi) + haem$hb_rap_m + log(haem$hb_pvr_calc) + log(haem$hb_pawp_m))

BIC(f1,f2,f3,f4,f5,f6,f7,f8)
anova(f1, f2)
anova(f1, f3)
anova(f1,f4)

#so adding in RAP gives most significance
#see if adding in extra helps
anova(f4, f8)

summary(f4)
resid_panel(f4,  plots = c("resid", "qq", "ls", "cookd"), smoother = TRUE)

#=============================================================================================
#run a correlation matrix for outcome variables
mat <- clin1[,c(2, 11, 24, 37, 150, 199, 203, 205, 206, 212, 313, 315)]

library(psych)

corr_mat=corr.test(mat, use='pairwise.complete.obs', method="s")
library(corrplot)

corrplot(corr_mat$r, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black", tl.cex = 0.7, tl.offset = 1, number.cex=0.75,
         tl.col = "black", insig = "pch",
         p.mat = corr_mat$p, sig.level = 0.05,
         title = "Spearman correlations of outcome variables")  

colnames(bic_scores_models_with_cci) <- bic_scores_models_with_cci[1,]
bic_scores_models_with_cci <- bic_scores_models_with_cci[-1,]

bc <- t(bic_scores_models_with_cci)
colnames(bc) <- bc[1,]
bc <- bc[-1,]
bc <- as.data.frame(bc)
bc[,1] <- as.numeric(bc[,1])
bc[,2] <- as.numeric(bc[,2])
bc[,3] <- as.numeric(bc[,3])
bc[,4] <- as.numeric(bc[,4])
bc[,5] <- as.numeric(bc[,5])
bc[,6] <- as.numeric(bc[,6])

bc <- as.data.frame(t(bc))

min1 = min(bc$PAWP)
min2 = min(bc$RAP)
min3 = min(bc$PVR)
min4 = min(bc$eGFR)
min5 = min(bc$FVC)
min6 = min(bc$`6MWD`)

bc[,1] <- bc[,1] - min1
bc[,2] <- bc[,2] - min2
bc[,3] <- bc[,3] - min3
bc[,4] <- bc[,4] - min4
bc[,5] <- bc[,5] - min5
bc[,6] <- bc[,6] - min6

write.csv2(bc, file='delta_BIC.csv')


