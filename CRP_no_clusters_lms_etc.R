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
load("C:/Users/Gebruiker/Downloads/data_clean (3).RData")
crp_df <- data.clean.cohort %>% dplyr::select(id, visit, id_cohort, cbt_inflammation_crp_mgpl, cbt_inflammation_scrp_mgpl)
crp_df <- as.data.frame(crp_df)
definitive.clusterings.3k.PAM.labels.3March2021 <- read.csv("~/Not PhD/My publications/Data new analyses March 2021/definitive clusterings 3k PAM labels 3March2021.csv", sep=";")
#===================================================================================
#Mortality script
#===================================================================================
#import the data
centre <- data.clean.cohort %>% dplyr::select('id', 'centre')
centre <- centre %>% filter(centre %in% c('Glasgow', 'Sheffield', 'Great Ormond Street', 'Lincoln', 'Papworth', 'Royal Brompton', 'Royal United Hospital Bath', 'Imperial and Hammersmith', 'Newcastle Freeman', 'Royal Free'))

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")
above_18 <- v4_clean_clinical_data_first_visit_15Dec23 %>% filter(age_diagnosis >=18)

crp_df <- crp_df %>% filter(id %in% above_18$id)


mort_df <-  above_18 %>% dplyr::select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'bs_bmi', 'DOB', 'sub_cause', 'sub_date')
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
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)
#=========================================================================================
#also import previous cluster labels
clusters_scaled_crp_time_pam_2k <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/clusters_scaled_crp_time_pam_2k.rds")
clusters_unscaled_crp_time_pam_2k <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/clusters_unscaled_crp_time_pam_2k.rds")


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

m_hl <- hl_crp %>% filter(visit <2)
write_rds(m_hl, 'crp_high_low_first_or_diagnostic_visit_v131Jan24.rds')

#now run mortality analysis
mdf <- merge(mort_df, m_hl)
mdf$above_5 <- as.factor(mdf$above_5)
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + bs_bmi + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)


#also plot without confidence intervals
ggsurvplot(sfit, data=mdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences based on first recorded CRP in PAH', xlim=c(0,15), legend.title='CRP >=5', break.x.by=2.5, legend.labs=c('First recorded CRP <=5', 'First recorded CRP >5'), palette = c('darkblue', 'darkred'))
 

#get clinical differences between groups
b <- merge(mdf, v4_clean_clinical_data_first_visit_15Dec23)
#so data for 808 pts

library(Publish)
a1 <- univariateTable(crp_bin ~ diagnosis_verified + sex + age_diagnosis + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + bs_bmi + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc + functional_class, Q.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b)
a2 <- summary(a1)

#now also do a survival analysis for just BMI
bmi_df <- mort_df %>% filter(age_diagnosis >=18)
bmi_df$weight <- ifelse(bmi_df$bs_bmi <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=18.5 & bmi_df$bs_bmi <25, 'healthy', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=25 & bmi_df$bs_bmi <30, 'overweight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=30 & bmi_df$bs_bmi <40, 'obese', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=40, 'severely_obese', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)

mdf <- bmi_df
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ weight, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ weight + sex + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)


ggsurvplot(sfit, data=mdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = FALSE, title='survival differences based on BMI', xlim=c(0,15), legend.title='BMI group', break.x.by=2.5, legend.labs=c('Healthy', 'Obese', 'Overweight', 'Severely obese', 'Underweight'), palette = c('darkgreen', 'darkorange', 'yellow', 'darkred', 'darkblue'))


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
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + bs_bmi + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)


#also plot without confidence intervals
ggsurvplot(sfit, data=mdf, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences based on CRP, after removal of outliers (*4 IQR from median)', xlim=c(0,15), legend.title='CRP >=5', break.x.by=2.5, legend.labs=c('First recorded CRP <=5', 'First recorded CRP >5'), palette = c('darkblue', 'darkred'))


#============================================================================
#now plot CRP trajectory if first CRP is >=5 vs below! 
hl_lab <- hl_crp %>% dplyr::select(id, above_5)


df_c <- data.clean.cohort %>% dplyr::select(id, visit, id_cohort, cbt_inflammation_crp_mgpl, cbt_inflammation_scrp_mgpl)
df_c <- as.data.frame(df_c)

#there is also HS-CRP data --> if CRP is NA use that 
df_c$crp_diff <- df_c$cbt_inflammation_crp_mgpl - df_c$cbt_inflammation_scrp_mgpl
summary(df_c$crp_diff)
#Differences are minimal
#I am aware of bias risk but roughly corresponds
df_c$cbt_inflammation_crp_mgpl <- ifelse(is.na(df_c$cbt_inflammation_crp_mgpl), df_c$cbt_inflammation_scrp_mgpl, df_c$cbt_inflammation_crp_mgpl)
df_c <- df_c %>% filter(!is.na(cbt_inflammation_crp_mgpl))

crp_dat <- df_c
crp_dat <- merge(crp_dat, hl_lab)
crp_dat$above_5 <- as.factor(crp_dat$above_5)

basic_plot <- crp_dat %>% group_by(above_5, visit) %>% dplyr::summarise(mean = mean(cbt_inflammation_crp_mgpl), stdev = sd(cbt_inflammation_crp_mgpl), meds = median(cbt_inflammation_crp_mgpl), iqr = 0.5*IQR(cbt_inflammation_crp_mgpl))
basic_plot$visit <- as.numeric(basic_plot$visit)
#only take up to visit 5
basic_plot <- basic_plot %>% filter(visit <6)


ggplot(data=basic_plot, aes(x=visit, y=meds, colour=as.factor(above_5))) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=meds-iqr, ymax=meds+iqr), width=.2) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  ggtitle('CRP levels in the cohort over time (median +- IQR)') +
  ylab('Median CRP level (mg/l)') +
  xlab('Visit') +
  scale_color_manual(values=c("darkblue", "darkred")) +
  guides(color = guide_legend(title = "First recorded CRP >=5"))


#plot geometric means too
library(plotrix)
gm_crp <- crp_dat %>% dplyr::group_by(above_5, visit) %>% dplyr::summarise(geo_mean = exp(mean(log(cbt_inflammation_crp_mgpl))), ci = plotrix::std.error(cbt_inflammation_crp_mgpl))
gm_crp <- gm_crp %>% filter(visit <=5)

plot_gm <- ggplot(gm_crp, aes(x=visit, y=geo_mean, colour=as.factor(above_5))) +
  geom_point() +
  geom_line() +
  xlab("Study visit") +
  ylab("CRP level (mg/ml)") +
  geom_errorbar(aes(ymin=geo_mean-ci, ymax=geo_mean+ci), width=.2) +
  ggtitle("Geometric means for CRP level, with 95% CI per group - UK cohort") + 
  scale_color_manual(values=c("darkblue", "darkred")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 20)) +
  guides(color = guide_legend(title = "CRP >=5"))


#=========================================
#now get repeated measures ANOVA for CRP

hist(crp_dat$cbt_inflammation_crp_mgpl)
hist(log(crp_dat$cbt_inflammation_crp_mgpl))

crp_anov <- crp_dat %>% filter(visit <6)
crp_anov$cbt_inflammation_crp_mgpl <- log(crp_anov$cbt_inflammation_crp_mgpl)
crp_anov$visit <- as.factor(crp_anov$visit)
crp_anov <- merge(crp_anov, hl_lab)


crp_anov_high <- crp_anov %>% filter(above_5 == 'yes')
crp_anov_low <- crp_anov %>% filter(above_5 == 'no')

library(datarium)
library(rstatix)
res.aov <- anova_test(data = crp_anov_high, dv = cbt_inflammation_crp_mgpl, wid = id, within = visit)
get_anova_table(res.aov)

res.aov <- anova_test(data = crp_anov_low, dv = cbt_inflammation_crp_mgpl, wid = id, within = visit)
get_anova_table(res.aov)

#repeat ANOVA with first visit removed
crp1 <- crp_anov
crp1$visit <- as.numeric(crp1$visit)
crp1 <- crp1 %>% group_by(id) %>% filter(visit !=min(visit))
crp1$id <- as.factor(crp1$id)
crp1$above_5 <- as.factor(crp1$above_5)
crp1$visit <- as.numeric(crp1$visit)

#get 2 way repeated measure time series
crp1 <- crp1 %>% filter(cbt_inflammation_crp_mgpl != -Inf)

#this doesn't work --> try something else:
model.aov <- aov(crp1$cbt_inflammation_crp_mgpl ~ 
                   crp1$above_5 * crp1$visit + 
                   Error(crp1$id/(crp1$above_5*crp1$visit)))

get_anova_table(res.aov)

#get clinical differences
b <- merge(m_hl[,c(1,7)], v4_clean_clinical_data_first_visit_15Dec23)
b$id <- as.factor(b$id)
b$pvri <- (b$hb_pap_m - b$hb_pawp_m)/b$ci

clin123 <- b

library(Publish)
a1 <- univariateTable(above_5 ~ diagnosis_verified + sex + age_diagnosis + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + bs_bmi + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + pvri + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc + functional_class, Q.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b)
a2 <- summary(a1)
b$above_5 <- as.factor(b$above_5)
wilcox.test(b$pvri ~ b$above_5)

b <- b %>% filter(ep_1_type_6mwt == 'corridor')
a1 <- univariateTable(above_5 ~ diagnosis_verified + sex + age_diagnosis + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + bs_bmi + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc + functional_class, Q.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b)
a2 <- summary(a1)

#make BMI factorial
hl_mort <- merge(mort_df, hl_lab)
bmi_df <- hl_mort
bmi_df$weight <- ifelse(bmi_df$bs_bmi <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=18.5 & bmi_df$bs_bmi <25, 'healthy', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=25 & bmi_df$bs_bmi <30, 'overweight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=30 & bmi_df$bs_bmi <40, 'obese', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=40, 'severely_obese', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)
hl_mort <- bmi_df
hl_mort$above_5 <- as.factor(hl_mort$above_5)

mdf <- hl_mort
mdf <- mdf %>% filter(!is.na(weight))
sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ above_5, data=mdf)
#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + weight + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)


#======================================
#EdB continued here on 17-01-2024
#compare comorbidity between groups
cleaned_comorbidity_dataV2_20Dec2023 <- readRDS("C:/Users/location/cleaned_comorbidity_dataV2_20Dec2023.rds")

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

clin4 <- merge(score_df2, clin123, by='id', all=T)
clin4$x <- ifelse(is.na(clin4$x), 0, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis <50, clin4$x + 0, clin4$x +1)
clin4$x <- ifelse(clin4$age_diagnosis >=60 & clin4$age_diagnosis <= 69, clin4$x + 1, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis >=70 & clin4$age_diagnosis <= 79, clin4$x + 2, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis >=80, clin4$x + 3, clin4$x)

clin4$above_5 <- as.factor(clin4$above_5)
clin4$CCI_adj_score <- clin4$x
clincci <- clin4 %>% select(id, CCI_adj_score)
saveRDS(clincci, 'CCI_score_per_patient.rds')

#===================================
#get complete CCI_score too
clin4 <- merge(score_df2, v4_clean_clinical_data_first_visit_15Dec23, by='id', all=T)
clin4$x <- ifelse(is.na(clin4$x), 0, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis <50, clin4$x + 0, clin4$x +1)
clin4$x <- ifelse(clin4$age_diagnosis >=60 & clin4$age_diagnosis <= 69, clin4$x + 1, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis >=70 & clin4$age_diagnosis <= 79, clin4$x + 2, clin4$x)
clin4$x <- ifelse(clin4$age_diagnosis >=80, clin4$x + 3, clin4$x)

clin4$CCI_adj_score <- clin4$x
clincci <- clin4 %>% select(id, CCI_adj_score)
saveRDS(clincci, 'CCI_score_per_patient_complete.rds')
####################################################################
clin4$above_5 <- ifelse(clin4$cbt_inflammation_crp_mgpl >5, 'yes', 'no')
clin4$above_5 <- ifelse(is.na(clin4$cbt_inflammation_crp_mgpl), NA, clin4$above_5)

clin4 <- clin4 %>% filter(!is.na(above_5))

u1 <- univariateTable(above_5 ~ x, data=clin4, compare.groups = TRUE, show.totals = TRUE, column.percent = TRUE)
u2 <- summary(u1)


p8 <- ggplot(clin4, aes(x=above_5, y=x, fill=above_5)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in Charlson Comorbidity Index between CRP groups') +
  ylab('CCI') +
  xlab('CRP group') +
  theme_bw() +
  stat_compare_means(method = 'wilcox.test') +
  scale_x_discrete(labels= labs)


#=============================
#add cox-PH with morbidity score
mdf <- merge(mort_df, m_hl)


#now also do a survival analysis for just BMI
bmi_df <- mort_df %>% filter(age_diagnosis >=18)
bmi_df$weight <- ifelse(bmi_df$bs_bmi <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=18.5 & bmi_df$bs_bmi <25, 'healthy', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=25 & bmi_df$bs_bmi <30, 'overweight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=30, 'obese_or_severly_obese', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)

mdf <- merge(bmi_df[,c(1,15)], mdf, by='id', all=T)
#mdf <- mdf %>% filter(diagnosis_verified== 'IPAH')
mdf <- merge(mdf, clin4[,c(1,11, 315, 203, 212)], by='id')
names(mdf)[[24]] <- 'CCI_score'
#also get tertiles
mdf$cci_tertile <- ifelse(mdf$CCI_score <0.1, 'Lowest Risk', 'Intermediate risk')
mdf$cci_tertile <- ifelse(mdf$CCI_score >3, 'Highest Risk', mdf$cci_tertile)
mdf$cci_tertile <- as.factor(mdf$cci_tertile)

#only retain patients who are adult and UK centres
mdf <- mdf %>% filter(!is.na(centre) & age_diagnosis.x >=18)
#736 retained
             
sobj <- Surv(mdf$surv_time, mdf$event, type='right')

#now transform data for scaling
mdf$age_scale <- scale(mdf$age_diagnosis.x)
mdf$pvr <- mdf$hb_pvr_calc/80
mdf$pvr_scale <- scale(mdf$pvr)
mdf$cci_scale <- scale(mdf$CCI_score)
mdf$rap_scale <- scale(mdf$hb_rap_m)
mdf$bmi_scale <- scale(mdf$bs_bmi)

mdf$above_5 <- as.factor(mdf$above_5)
mdf$sex <- as.factor(mdf$sex)

#do cox
cox <- coxph(sobj ~ above_5 + sex + age_scale + bmi_scale + cci_scale + pvr_scale + rap_scale, data=mdf)
ggforest(cox, data=mdf)


#save now
sc <- summary(cox)
sc2 <- as.data.frame(sc$coefficients)
saveRDS(sc2,file='coefficients_SCALED_cox_per_CRP_cohort.rds')


#based on CW's question --> add in smoking
query_smoking_dump_10_7_23 <- read.delim("~/PhD/Projects/CRP ~ survival and BMI/query_smoking_dump_10_7_23.txt")
smoke <- query_smoking_dump_10_7_23
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_history == 'NULL', 'no', smoke$cfh_smoking_history)
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_history == 'UNK', NA, smoke$cfh_smoking_history)
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_current == 'yes', 'Current smoker', smoke$cfh_smoking_history)
smoke$cfh_smoking_history <- ifelse(smoke$cfh_smoking_history == 'yes', 'Smoking history', smoke$cfh_smoking_history)

smoke$`Smoking history` <- as.factor(smoke$cfh_smoking_history)

pvr_smoke <- merge(smoke[,c(2,6)], v4_clean_clinical_data_first_visit_15Dec23, by.y='id', by.x='study_subject_oid')
pvr_smoke <- pvr_smoke %>% filter(!is.na(`Smoking history`))

sample_size <- function(x) {
  return(c(y = max(x) + 0.5, label = length(x)))
}

pvr_smoke$smoke <- ifelse(pvr_smoke$`Smoking history` %in% c('Current smoker', 'Smoking history'), 'Smoker', 'Non smoker')


ggplot(pvr_smoke, aes(x=smoke, y=hb_pvr_calc, fill=smoke)) +
         geom_violin() +
         geom_boxplot(width=0.1) +
         stat_compare_means(method='wilcox.test') +
  ggtitle('Effect of smoking on PVR') +
  xlab('Smoking history') +
  #stat_summary(fun.data = sample_size, geom = "text", vjust = 40, color = "black") +
  ylab('PVR (dyne/sec/cm^-5)')

pvr_smoke2 <- pvr_smoke %>% filter(`Smoking history` %in% c('Current smoker', 'Smoking history'))

ggplot(pvr_smoke2, aes(x=`Smoking history`, y=hb_pvr_calc, fill=`Smoking history`)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  stat_compare_means(method='wilcox.test') +
  ggtitle('Effect of smoking on PVR') +
  xlab('Smoking history') +
  #stat_summary(fun.data = sample_size, geom = "text", vjust = 40, color = "black") +
  ylab('PVR (dyne/sec/cm^-5)')


#merge with CoxPH data
mdf1234 <- merge(mdf, smoke[,c(2,6)], by.x='id', by.y='study_subject_oid')

sobj <- Surv(mdf1234$surv_time, mdf1234$event, type='right')

levels(mdf1234$`Smoking history`)
mdf1234$`Smoking history` <- relevel(mdf1234$`Smoking history`, ref=3)

cox <- coxph(sobj ~  above_5 + sex + age_scale + bmi_scale + cci_scale + pvr_scale + rap_scale + `Smoking history`, data=mdf1234)
ggforest(cox, data=mdf1234)





#now run cox-pcox#now run cox-ph for age + sex + bmi
cox <- coxph(sobj ~ above_5 + sex + weight + cci_tertile, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ above_5 + sex + bs_bmi + cci_tertile, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi + age_diagnosis, data=mdf)
ggforest(cox, data=mdf)


cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi + age_diagnosis + ci, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi + age_diagnosis + hb_rap_m + hb_pvr_calc + hb_pawp_m, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi + age_diagnosis + hb_pvr_calc, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi + hb_pap_m + hb_cardiac_output_value_1 + hb_rap_m + hb_pawp_m, data=mdf)
ggforest(cox, data=mdf)


cox <- coxph(sobj ~ above_5 + sex + CCI_score + bs_bmi + hb_rap_m + hb_pvr_calc, data=mdf)
ggforest(cox, data=mdf)

#include smoking status
query_smoking_dump_10_7_23 <- read.delim("~/PhD/Projects/CRP ~ survival and BMI/query_smoking_dump_10_7_23.txt")
smoke <- merge(query_smoking_dump_10_7_23, crp_high_low_first_or_diagnostic_visit_v131Jan24, by.x='study_subject_oid', by.y='id')
smoke2 <- smoke %>% select(study_subject_oid, cfh_smoking_history, cfh_smoking_current, above_5)
smoke2 <- unique(smoke2)
smoke2$cfh_smoking_history <- ifelse(smoke2$cfh_smoking_history == 'NULL', 'no', smoke2$cfh_smoking_history)
smoke2$cfh_smoking_history <- ifelse(smoke2$cfh_smoking_history == 'UNK', NA, smoke2$cfh_smoking_history)
smoke2$cfh_smoking_history <- as.factor(smoke2$cfh_smoking_history)
chisq.test(smoke2$above_5, smoke2$cfh_smoking_history)

smoke2$smoke_now <- ifelse(smoke2$cfh_smoking_history == 'yes', 'past_smoker', 'non_smoker')
smoke2$smoke_now <- ifelse(smoke2$cfh_smoking_current == 'yes', 'current_smoker', smoke2$smoke_now)
smoke2$smoke_now <- as.factor(smoke2$smoke_now)
chisq.test(smoke2$above_5, smoke2$smoke_now)

mdf3 <- merge(mdf, smoke2[,c(1,5)], by.x='id', by.y='study_subject_oid')
s1 <- Surv(mdf3$surv_time, mdf3$event, type='right')

cox <- coxph(s1 ~ above_5 + sex + smoke_now + CCI_score + bs_bmi + age_diagnosis + hb_rap_m, data=mdf3)
ggforest(cox, data=mdf3)

#now with just smoking and CCI
cox <- coxph(s1~ smoke_now + CCI_score,data=mdf3)
ggforest(cox, data=mdf3)

cox <- coxph(s1~ smoke_now,data=mdf3)
ggforest(cox, data=mdf3)
#=====================================================


cox <- coxph(sobj2 ~ above_5 + definitive_labels, data=mdf2)
ggforest(cox, data=mdf2)

summary(aov(cbt_inflammation_crp_mgpl ~ definitive_labels, data=mdf2))
#==================================================
#plot clinical differences
labs <- c('<=5', '>5')
p1 <- ggplot(clin123, aes(x=above_5, y=hb_pawp_m, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('PAWP') +
  ylab('Wedge pressure (mmHg)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

p2 <- ggplot(clin123, aes(x=above_5, y=hb_rap_m, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('RAP') +
  ylab('RAP (mmHg)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

p3 <- ggplot(clin123, aes(x=above_5, y=bs_bmi, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('BMI') +
  ylab('BMI (kg/m^2)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

p4 <- ggplot(clin123, aes(x=above_5, y=egfr_mdrd, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('eGFR') +
  ylab('eGFR (MDRD;ml/min/1.73m^2)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

b <- clin123 %>% filter(ep_1_type_6mwt == 'corridor')

p5 <- ggplot(b, aes(x=above_5, y=ep_1_distance_meters, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('6MWD') +
  ylab('6MWD (m)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)


p6 <- ggplot(clin123, aes(x=above_5, y=lf_fvc_pc, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('FVC') +
  ylab('FVC (%pred)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)


ggplot(clin123, aes(x=above_5, y=cfe_rest_spo2, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in resting SpO2 between CRP groups') +
  ylab('SpO2 (%)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

clin123$functional_class <- as.factor(clin123$functional_class)
fc <- clin123 %>% dplyr::group_by(above_5, functional_class) %>% dplyr::summarise(n())
fc <- fc %>% filter(!is.na(functional_class))
names(fc)[[3]] <- 'number'
fc$percentage <- ifelse(fc$above_5 == 'yes', (fc$number/381)*100, (fc$number/668)*100)


ggplot(fc, aes(x=functional_class, y=percentage, fill=functional_class)) + 
  geom_bar(stat = 'identity') +
  ggtitle('Difference in FC between CRP groups') +
  ylab('Percentage of group') +
  xlab('WHO functional class') +
  facet_wrap(~above_5)

p7 <- ggplot(clin123, aes(x=above_5, y=hb_pvr_calc, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('PVR') +
  ylab('PVR (dynes/sec)') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

ggplot(clin123, aes(x=above_5, y=pvri, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in PVRi between CRP groups') +
  ylab('Indexed PVR') +
  xlab('CRP group') +
  theme_bw() +
  scale_x_discrete(labels= labs)

library(patchwork)
p1 + p2 + p3 +p4+p5+p6+p7 + plot_layout(ncol=7)

#import the AABs
library(readxl)
log.transformed.autoAbs.of.ctrls.and.pts.v1.3March2021 <- read.csv2("~/Not PhD/My publications/Data new analyses March 2021/log-transformed autoAbs of ctrls and pts v1 3March2021.csv")
autoimmunity_dataframe <- read_excel("~/Not PhD/TBR PAH cohort study/Origional data with minor transformations/Original data/autoimmunity.dataframe.xlsx")
autoimmunity_dataframe <- as.data.frame(autoimmunity_dataframe)
colnames(autoimmunity_dataframe) <- autoimmunity_dataframe[1,]
autoimmunity_dataframe <- autoimmunity_dataframe[-1,]
ai <- autoimmunity_dataframe %>% select(id_oc, id_cohort.x)
abs <- merge(log.transformed.autoAbs.of.ctrls.and.pts.v1.3March2021, ai, by.x='sample', by.y='id_cohort.x')

ab2 <- merge(abs, clin123[,c(1,2)], by.x = 'id_oc', by.y='id')
ab2 <- unique(ab2)

names(ab2) <- c('id', 'sample', 'x', 'Cardiolipin', 'CENP-B', 'H2a/F2a & H4a/F2a1', 'Histone IIA', 'Jo-1', 'La/SS-B', 'Mi-2b', 'MPO', 'Proteinase-3', 'PDH', 'RNP complex', 'Ro/SS-A', 'Scl-34', 'Scl-70', 'Smith', 'Thyroglobulin', 'TPO', 'Transglutaminase', 'U1-snRNP68', 'sex', 's', 'age', 'pah', 'above_5')
#now group AABs
ab2 <- ab2 %>% pivot_longer(4:22, names_to='Autoantibody', values_to = 'AAB_level')

library(ggpubr)

ggplot(ab2, aes(x = above_5, y = AAB_level, fill = above_5)) + 
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 4) +
  ggtitle('Difference in Autoantibody level between CRP groups') +
  ylab('Log-normalised AAB level') +
  xlab('CRP >5') +
  theme_bw() +
  #stat_compare_means(method = 'wilcoxon.test', label = "p.format", vjust = 1) +
  facet_wrap(~ Autoantibody, scales = "free_y") +
  scale_fill_manual(values = c('darkblue', 'darkred'))

results <- ab2 %>%
  group_by(Autoantibody) %>%
  do(tidy(wilcox.test(AAB_level ~ above_5, data = .)))

# Print the results
print(results)




#=================================================================================
#now repeat baseline table!
library(Publish)
library(broom)
above_18_t <- merge(above_18, hl_lab, by='id')
#968 patients


a1 <- univariateTable(above_5 ~ diagnosis_verified + sex + age_diagnosis + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + bs_bmi + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc + functional_class + cbt_card_ntprobnp_ngpl, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = above_18_t)
a2 <- summary(a1)

b <- above_18_t

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


group_var = "above_5"  
output <- lapply(myvar_ch, function(x) do_chisq(b, x, group_var))
combined_output1 <- bind_rows(output)


chisq.test(b$diagnosis_verified, b$above_5)
#0.1368
chisq.test(b$hv_vasodilator_responder, b$above_5)
#0.0433

combined_output1$p.value <- ifelse(combined_output1$var == 'diagnosis_verified', 0.1368, combined_output1$p.value)
combined_output1$p.value <- ifelse(combined_output1$var == 'hv_vasodilator_responder', 0.0433, combined_output1$p.value)


#make function for mann-whithney-U test

do_mw <- function(myvar) {
  mw_output <- wilcox.test(b[[myvar]] ~ b$above_5)
  mw_output <- (tidy(mw_output))
  mw_output$var <- myvar
  return(mw_output)
}

do_mw('age_diagnosis')

#now do for all
myvar_mw <- c('age_diagnosis', 'cbt_thyr_tsh_mupl', 'cbt_thyr_freet4_pmolpl', 'bs_bmi', 'egfr_mdrd', 'cbt_haem_platelets_x10e9pl', 'cbt_haem_hb_gpl', 'cfe_rest_spo2', 'cfe_heart_rate', 'hb_pawp_m', 'hb_pap_d',  'hb_pap_m', 'hb_pvr_calc', 'hb_rap_m', 'hb_cardiac_output_value_1', 'hb_cardiac_index_value_1', 'lf_fev1_pc', 'lf_fvc_pc', 'lf_kco_pc', 'cbt_card_ntprobnp_ngpl')
myvar = myvar_mw
output <- lapply(myvar, do_mw)
combined_output2 <- bind_rows(output)

pvals <- rbind(combined_output1[,c(1,2,4,5)], combined_output2[,c(1:3,5)])

pvals$fdr <- p.adjust(pvals$p.value, method='fdr')
pvals$fdr_round <- round(pvals$fdr, digits=10)

#quick repeat with just corridor walk
b2 <- b %>% filter(ep_1_type_6mwt == 'corridor')
c1 <- univariateTable(above_5 ~ ep_1_distance_meters, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b2)
c2 <- summary(c1)

wilcox.test(b2$ep_1_distance_meters ~ b2$above_5)
#pval = 1.224e-08

pvals[25,] <- NA
pvals$p.value <- ifelse(is.na(pvals$var), 1.224e-08, pvals$p.value)
pvals$var <- ifelse(is.na(pvals$var), 'ep_1_distance_meters', pvals$var)

#get pval for CCI scores (from earlier plot)
pvals[26,] <- NA
pvals$p.value <- ifelse(is.na(pvals$var), 0.017, pvals$p.value)
pvals$var <- ifelse(is.na(pvals$var), 'CCI_score', pvals$var)

pvals$fdr <- p.adjust(pvals$p.value, method='fdr')
pvals$fdr_round <- round(pvals$fdr, digits=10)


#get a CCI table
CCI_score_per_patient_complete <- readRDS("C:/Users/location/CCI_score_per_patient_complete.rds")
b3 <- merge(b, CCI_score_per_patient_complete, by='id')

d1 <- univariateTable(above_5 ~ CCI_adj_score, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b3)
d2 <- summary(d1)

ac <- rbind(a2, d2)


#p8 <- ggplot(b3, aes(x=above_5, y=CCI_adj_score, fill=above_5)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in Charlson Comorbidity Index between CRP groups') +
  ylab('CCI') +
  xlab('CRP group') +
  theme_bw() +
  stat_compare_means(method = 'wilcox.test') +
  scale_x_discrete(labels= labs)

#now merge these data with the univariate table

pv2 <- pvals %>% select(var, p.value, fdr_round)

ac$increasing_numbers <- seq(1, nrow(ac))


a3 <- merge(ac, pv2, by.x='Variable', by.y='var', all=T)

a3 <- a3[order(a3$increasing_numbers), ]

write.csv2(a3, file='univar_table_baseline_crp_high_low.csv')
