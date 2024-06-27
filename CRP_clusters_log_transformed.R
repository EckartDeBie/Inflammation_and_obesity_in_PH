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
load("C:/filepath.RData")
#==========================================================================
#now do survival analysis
#import the data
centre <- data.clean.cohort %>% select('id', 'centre')
centre <- centre %>% filter(centre %in% c('Glasgow', 'Sheffield', 'Great Ormond Street', 'Lincoln', 'Papworth', 'Royal Brompton', 'Royal United Hospital Bath', 'Imperial and Hammersmith', 'Newcastle Freeman', 'Royal Free'))

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/filepath.rds")

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

library(survival)
library(survminer)

#=======================================================================
#Now get CRP trajectories in nice ggplot with ranges
#=======================================================================
over18 <- v4_clean_clinical_data_first_visit_15Dec23 %>% filter(age_diagnosis >=18)
data.clean.cohort <- data.clean.cohort %>% filter(id %in% over18$id)
df_c <- data.clean.cohort %>% select(id, visit, id_cohort, cbt_inflammation_crp_mgpl, cbt_inflammation_scrp_mgpl)
df_c <- as.data.frame(df_c)

#there is also HS-CRP data --> if CRP is NA use that 
df_c$crp_diff <- df_c$cbt_inflammation_crp_mgpl - df_c$cbt_inflammation_scrp_mgpl
summary(df_c$crp_diff)
#Differences are minimal
#I am aware of bias risk but roughly corresponds
df_c$cbt_inflammation_crp_mgpl <- ifelse(is.na(df_c$cbt_inflammation_crp_mgpl), df_c$cbt_inflammation_scrp_mgpl, df_c$cbt_inflammation_crp_mgpl)
df_c <- df_c %>% filter(!is.na(cbt_inflammation_crp_mgpl))

crp_dat <- df_c
crp_for_later <- crp_dat

#run some DTW time series clustering on the data

#import the packages
library(readxl)
library(data.table)
library(magrittr)
library(tidyverse)
#probably required later on
library(ConsensusClusterPlus)
library(pheatmap)
#library(corrplot)
library(factoextra)
library(dtw)
library(dtwclust)

#used background information on DTW from: https://rpubs.com/esobolewska/dtw-time-series
#general approach
#1) take dtw with dtw distance, and use PAM & Km
#2) investigate optimal K and algorithm

#get the data in shape
time_df <- df_c[,c(1,2,4)] %>% pivot_wider(names_from = visit, values_from = cbt_inflammation_crp_mgpl)
#given missingness only include first 5 visits
time_df <- time_df[,c(1:6)]
time_df <- unique(time_df)

time_df <- as.data.frame(time_df)
rownames(time_df) <- time_df$id
time_df <- time_df[,-c(1)]
t2 <- time_df
#the data is in the correct format

#DTW clust does not accept missing values! 
#means impute
#if visit 4 has nothing recorded --> take visit 3's value
time_df$'4' <- ifelse(is.na(time_df$'4'), time_df$'3', time_df$'4')
#if visit 0 has no CRP --> take visit 1's CRP
time_df$'0' <- ifelse(is.na(time_df$'0'), time_df$'1', time_df$'0')
time_df$'3' <- ifelse(is.na(time_df$'3'), (time_df$'2' + time_df$'4')/2, time_df$'3')
time_df$'2' <- ifelse(is.na(time_df$'2'), (time_df$'1' + time_df$'3')/2, time_df$'2')
time_df$'1' <- ifelse(is.na(time_df$'1'), (time_df$'0' + time_df$'2')/2, time_df$'1')
#check how many pts are kept in with all other NAs removed
try <- na.omit(time_df)
#178/1230 values imputed = 14.5%
t3 <- try
t3$tname <- 'x'
t3$id <- rownames(t3)
t2$id <- rownames(t2)
t2 <- merge(t2, t3[,6:7], by='id')

lapply(t2, summary)
23 + 22 + 22 + 10 + 87
5*236
#164/1180 values imputed = 13.9%

try2 <- try
#236 patients are kept in! 
#now make sure the data is log-transformed
try$id <- rownames(try)
try[,1:5] <- lapply(try[,1:5], log)
try <- as.data.frame(try[,-6])


#===================================================================================
#Do sense checks

pam_t1 <- tsclust(try, type="partitional", k=2L:8L, distance="dtw", centroid="pam")

sense_checks <- lapply(pam_t1, cvi)
#for this, the following applies:
#Sil, SF, CH, D need to be max, and CH, DBstar, COP need to be minimal for optimal clustering
#the CVI function is an easy way to call all FVIZ_nbclust and others with minimal computation

#make sure to get a long list and plot everything
list_checks <- data.table::rbindlist(lapply(sense_checks, as.data.table, keep.rownames=TRUE), use.names = TRUE, idcol = TRUE)

#now find optimal clustering
ggplot(list_checks, aes(x=.id, y=V2)) +
  geom_bar(stat='identity') +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7)) +
  facet_wrap(~V1, scales = 'free')


#repeat this for HC
hc_t1 <- tsclust(try, type = "hierarchical", k = 2L:8L, distance = "dtw")
sense_checks <- lapply(hc_t1, cvi)
#for this, the following applies:
#Sil, SF, CH, D need to be max, and CH, DBstar, COP need to be minimal for optimal clustering

#the CVI function is an easy way to call all FVIZ_nbclust and others with minimal computation

#make sure to get a long list and plot everything
list_checks <- data.table::rbindlist(lapply(sense_checks, as.data.table, keep.rownames=TRUE), use.names = TRUE, idcol = TRUE)

#plot
ggplot(list_checks, aes(x=.id, y=V2)) +
  geom_bar(stat='identity') +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7)) +
  facet_wrap(~V1, scales = 'free')

#========================================================
#Do actual clustering
#continue with PAM
#now cluster with DTW for this!
def_pam_t1 <- tsclust(try, type="partitional", seed = 123, k=2L, distance="dtw", clustering="pam")

#plot the clusters with each group member included
plot(def_pam_t1, type = "sc")

#==============================================
#now cluster with DTW for this!
hc_t1 <- tsclust(try, type="hierarchical", seed = 123, k=2L, distance="dtw", clustering="hierarchical")

#plot the clusters with each group member included
plot(hc_t1, type = "sc")
#=======================================

#see to which cluster each patient belongs
clusters_def <- t(cbind(try2[,0], cluster = def_pam_t1@cluster))
clusters_def <- as.data.frame(clusters_def)
clusters_def <- pivot_longer(clusters_def, 1:236)

try2$name <- rownames(try2)
df3 <- merge(try2, clusters_def)
df3 <- df3 %>% pivot_longer(2:6, names_to = 'visit', values_to = 'cbt_inflammation_crp_mgpl')

basic_plot <- df3 %>% group_by(visit, value) %>% summarise(mean = mean(cbt_inflammation_crp_mgpl), stdev = sd(cbt_inflammation_crp_mgpl), meds = median(cbt_inflammation_crp_mgpl), iqr = 0.5*IQR(cbt_inflammation_crp_mgpl))

df4 <- clusters_def

saveRDS(unique(df4), file='log_transformed_crp_levels_clusters.rds')

ggplot(data=basic_plot, aes(x=as.numeric(visit), y=meds, colour=as.factor(value))) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=meds-iqr, ymax=meds+iqr), width=.2) +
  scale_x_continuous(breaks = c(0,1,2,3,4)) +
  ggtitle('CRP levels in the cohort over time (median + IQR)') +
  ylab('Median CRP level (mg/l)') +
  xlab('Visit') +
  scale_color_manual(values=c("darkred", "darkblue")) +
  guides(color = guide_legend(title = "DTW cluster after log-scaling CRP"))


#get geometric means
library(plotrix)
df3$value <- as.factor(df3$value)
df3$visit <- as.numeric(df3$visit)
gm_crp <- df3 %>% dplyr::group_by(value, visit) %>% dplyr::summarise(geo_mean = exp(mean(log(cbt_inflammation_crp_mgpl))), ci = plotrix::std.error(cbt_inflammation_crp_mgpl))

gm_crp$value <- ifelse(gm_crp$value == '1', 'Normal CRP cluster', 'High CRP cluster')
gm_crp$value <- as.factor(gm_crp$value)

plot_gm <- ggplot(gm_crp, aes(x=visit, y=geo_mean, colour=as.factor(value))) +
  geom_point() +
  geom_line() +
  xlab("Study visit") +
  ylab("CRP level (mg/ml)") +
  geom_errorbar(aes(ymin=geo_mean-ci, ymax=geo_mean+ci), width=.2) +
  ggtitle("Geometric means with 95% CI for each CRP cluster") + 
  scale_color_manual(values=c("darkred", "darkblue")) +
  theme_bw() +
  guides(color = guide_legend(title = "CRP cluster"))


#=======================================================
#see if trajectory is stable
crp_anov <- df3

crp_anov_high <- crp_anov %>% filter(value == 1)
crp_anov_low <- crp_anov %>% filter(value == 2)

library(datarium)
library(rstatix)
res.aov <- anova_test(data = crp_anov_high, dv = cbt_inflammation_crp_mgpl, wid = name, within = visit)
get_anova_table(res.aov)

res.aov <- anova_test(data = crp_anov_low, dv = cbt_inflammation_crp_mgpl, wid = name, within = visit)
get_anova_table(res.aov)

crp1 <- crp_anov
crp1$visit <- as.numeric(crp1$visit)
crp1 <- crp1 %>% group_by(name) %>% filter(visit !=min(visit))
crp1$name <- as.factor(crp1$name)
crp1$value <- as.factor(crp1$value)
crp1$visit <- as.numeric(crp1$visit)

#get 2 way repeated measure time series
#first log-transform
crp1$cbt_inflammation_crp_mgpl <- log(crp1$cbt_inflammation_crp_mgpl)
crp1 <- crp1 %>% filter(cbt_inflammation_crp_mgpl != -Inf)

#this doesn't work --> try something else:
model.aov <- aov(crp1$cbt_inflammation_crp_mgpl ~ 
                   crp1$value * crp1$visit + 
                   Error(crp1$name/(crp1$value*crp1$visit)))

summary(model.aov)
#==========================================================================================
#do mortality analysis
mort_df

lab_mort <- merge(unique(df3[,1:2]), mort_df, by.x = 'name', by.y='id')
#241 patients kept in for survival analysis
lab_mort$value <- as.factor(lab_mort$value)

sobj <- Surv(lab_mort$surv_time, lab_mort$event, type='right')
sfit <- survfit(sobj ~ value, data=lab_mort)

#now run cox-ph for age + sex
cox <- coxph(sobj ~ value + sex + age_diagnosis, data=lab_mort)
ggforest(cox, data=lab_mort)
#no effect

#also plot without confidence intervals
ggsurvplot(sfit, data=lab_mort, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival differences based on CRP', xlim=c(0,15), legend.title='CRP cluster', break.x.by=2.5, legend.labs=c('high CRP cluster', 'low CRP cluster'), palette = c('darkred', 'darkblue'))

#import the clinical data
v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/filepath.rds")

b <- merge(lab_mort, v4_clean_clinical_data_first_visit_15Dec23)
#so data for 241 pts

#=============================================================================
#get a power calculation for this dataset
#code adapted from:https://shariq-mohammed.github.io/files/cbsa2019/2-power-and-sample-size.html
#this code doesn't work (due to immortal time bias in our dataset)......
library(Hmisc)
cpower(tref = 15, mc=38/133, accrual =0, r=36 , tmin = 15, noncomp.c = 0, noncomp.i = 0, alpha=0.05, nc=133, ni=108, pr=TRUE)

res.sum <- surv_summary(sfit, data = lab_mort)
summary(sfit, 15)
summary(sfit, 10)

#now cut at 10 years
cpower(tref = 10, mc= 26/133, accrual =0, r=36 , tmin = 10, noncomp.c = 0, noncomp.i = 0, alpha=0.05, nc=133, ni=108, pr=TRUE)


#Repeat at 5 yrs
summary(sfit, 5)
cpower(tref = 5, mc= 12/133, accrual =0, r=36 , tmin = 10, noncomp.c = 0, noncomp.i = 0, alpha=0.05, nc=133, ni=108, pr=TRUE)
#==============================================================================
load("C:/filepath.RData")
sel <- data.checked %>% select(id, cbt_card_bnp_ngpl)
b <- merge(b, sel, by='id')


library(Publish)
library(broom)
a1 <- univariateTable(value ~ diagnosis_verified + sex + age_diagnosis + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + bs_bmi + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc + functional_class + cbt_card_ntprobnp_ngpl + cbt_card_bnp_ngpl, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b)
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


group_var = "value"  
output <- lapply(myvar_ch, function(x) do_chisq(b, x, group_var))
combined_output1 <- bind_rows(output)


chisq.test(b$diagnosis_verified, b$value)

#make function for mann-whithney-U test

do_mw <- function(myvar) {
  mw_output <- wilcox.test(b[[myvar]] ~ b$value)
  mw_output <- (tidy(mw_output))
  mw_output$var <- myvar
  return(mw_output)
}

do_mw('age_diagnosis')

#now do for all
myvar_mw <- c('age_diagnosis', 'cbt_thyr_tsh_mupl', 'cbt_thyr_freet4_pmolpl', 'bs_bmi', 'egfr_mdrd', 'cbt_haem_platelets_x10e9pl', 'cbt_haem_hb_gpl', 'cfe_rest_spo2', 'cfe_heart_rate', 'hb_pawp_m', 'hb_pap_d',  'hb_pap_m', 'hb_pvr_calc', 'hb_rap_m', 'hb_cardiac_output_value_1', 'hb_cardiac_index_value_1', 'lf_fev1_pc', 'lf_fvc_pc', 'lf_kco_pc', 'cbt_card_ntprobnp_ngpl', 'cbt_card_bnp_ngpl')
myvar = myvar_mw
output <- lapply(myvar, do_mw)
combined_output2 <- bind_rows(output)

pvals <- rbind(combined_output1[,c(1,2,4,5)], combined_output2[,c(1:3,5)])

chisq.test(b$diagnosis_verified, b$value)
pvals$p.value <- ifelse(pvals$var == 'diagnosis_verified', 0.7052, pvals$p.value)

chisq.test(b$hv_vasodilator_responder, b$value)
pvals$p.value <- ifelse(pvals$var == 'hv_vasodilator_responder', 1, pvals$p.value)
pvals$fdr <- p.adjust(pvals$p.value, method='fdr')
pvals$fdr_round <- round(pvals$fdr, digits=10)

#quick repeat with just corridor walk
b2 <- b %>% filter(ep_1_type_6mwt == 'corridor')
c1 <- univariateTable(value ~ ep_1_distance_meters, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = b2)
c2 <- summary(c1)

wilcox.test(b2$ep_1_distance_meters ~ b2$value)
#pval = 0.0004523

pvals[26,] <- NA
pvals$p.value <- ifelse(is.na(pvals$var), 0.0004523, pvals$p.value)
pvals$var <- ifelse(is.na(pvals$var), 'ep_1_distance_meters', pvals$var)

pvals$fdr <- p.adjust(pvals$p.value, method='fdr')
pvals$fdr_round <- round(pvals$fdr, digits=10)


#now merge these data with the univariate table

pv2 <- pvals %>% select(var, p.value, fdr_round)

a2$increasing_numbers <- seq(1, nrow(a2))


a3 <- merge(a2, pv2, by.x='Variable', by.y='var', all=T)

a3 <- a3[order(a3$increasing_numbers), ]

write.csv2(a3, file='univar_table_log_crp_clusters.csv')




#check effect of cluster on AAB cluster
ids <- data.clean.cohort %>% select(id, id_cohort)
ids <- merge(ids, definitive.clusterings.3k.PAM.labels.3March2021, by.x='id_cohort', by.y='X')
ids <- merge(clusters_def, ids, by.x='name', by.y='id')
ids <- ids %>% select(name, value, definitive_labels)
ids <- unique(ids)
#171 patients retained
ids$value <- as.factor(ids$value)
ids$definitive_labels <- as.factor(ids$definitive_labels)
c1 <- univariateTable(definitive_labels ~ value, data=ids, column.percent = TRUE, show.totals = TRUE, compare.groups = TRUE)
c2 <- summary(c1)

#==================================================================
#see if clusters are comparable; adjusted Rand index --> import data
clusters_scaled_crp_time_pam_2k <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/clusters_scaled_crp_time_pam_2k.rds")
clusters_unscaled_crp_time_pam_2k <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/clusters_unscaled_crp_time_pam_2k.rds")

lclust <- df3 %>% select(name, value)
names(lclust) <- c('name', 'log_clust')

sclust <- clusters_scaled_crp_time_pam_2k
names(sclust) <- c('name', 'scaled_clust')

uclust <- clusters_unscaled_crp_time_pam_2k
names(uclust) <- c('name', 'unscaled_clust')

#merge everything into one df
cluster_rand <- merge(lclust, sclust, by='name')
cluster_rand <- merge(cluster_rand, uclust, by='name')
cluster_rand <- unique(cluster_rand)
#all observations have been retained

#get rand index
library(mclust)
adjustedRandIndex(cluster_rand$log_clust, cluster_rand$unscaled_clust)
adjustedRandIndex(cluster_rand$log_clust, cluster_rand$scaled_clust)
adjustedRandIndex(cluster_rand$unscaled_clust, cluster_rand$scaled_clust)

cluster_rand$log_is_scaled <- cluster_rand$log_clust == cluster_rand$scaled_clust
cluster_rand$log_is_unscaled <- cluster_rand$log_clust == cluster_rand$unscaled_clust
cluster_rand$scaled_is_unscaled <- cluster_rand$scaled_clust == cluster_rand$unscaled_clust

lapply(cluster_rand, summary)
