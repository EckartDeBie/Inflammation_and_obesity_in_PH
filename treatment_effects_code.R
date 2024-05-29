#assess treatment effects per groups
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
library(pROC)

#import the data
load("C:/Users/Gebruiker/Downloads/data_clean (3).RData")

#select relevant cols
treat_effect <- data.clean.cohort %>% select(id, visit, ep_1_type_6mwt, ep_1_distance_meters, functional_class)
#now get only visit 0 and 1 for treatment response

treat_effect <- treat_effect %>% filter(visit %in% c(0,1))
treat2 <- treat_effect %>% pivot_wider(names_from = visit, values_from = c(ep_1_distance_meters, ep_1_type_6mwt, functional_class))

#make sure the walk is the same type
treat2$walk_same <- treat2$ep_1_type_6mwt_0 == treat2$ep_1_type_6mwt_1
treat2$walk_diff_meters <- treat2$ep_1_distance_meters_1 - treat2$ep_1_distance_meters_0
treat2$walk_diff_meters <- ifelse(treat2$walk_same == 'FALSE', NA, treat2$walk_diff_meters)
treat2$fc_diff <- treat2$functional_class_1 - treat2$functional_class_0

#now look at distributions
hist(treat2$walk_diff_meters)
#nice normal distribution
hist(treat2$fc_diff)
#also nice normal distribution
#get % for vars
treat2$walk_diff_percent <- (treat2$walk_diff_meters/treat2$ep_1_distance_meters_0)*100
treat2$fc_diff_percent <- (treat2$fc_diff/treat2$functional_class_0)*100

#plot for distributions
hist(treat2$walk_diff_percent)
#obvious skew!
hist(treat2$fc_diff_percent)
#also a clear skew!

#import the bmi data
v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location.rds")
bmi <- v4_clean_clinical_data_first_visit_15Dec23 %>% select(id, bs_bmi)
treat3 <- merge(treat2, bmi, by='id')

saveRDS(treat3, file='treatment_effects_file_v1_March2024.rds')


BMI_category_per_patient <- readRDS("~/PhD/Projects/CRP ~ survival and BMI/BMI_category_per_patient.rds")
t2 <- merge(treat2, BMI_category_per_patient, by='id')
t2 <- t2 %>% filter(walk_same == 'TRUE')

#repeat for BMI groups
levels(t2$weight)

#relevel
desired_order <- c("obesity class-III", "obesity class-II","obesity class-I","pre-obesity","normal weight","underweight")

# Reorder the levels of the factor
t2$weight <- factor(t2$weight, levels = desired_order)
#only take 6MWD
t2 <- t2 %>% filter(ep_1_type_6mwt_0 == 'corridor')
#238 patients left!
t2 <- t2 %>% filter(!is.na(weight))
#230 left

ggplot(t2, aes(x=weight, y=walk_diff_meters)) + 
  geom_boxplot(aes(fill = weight), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in 6MWD improvement between BMI groups in PAH - UK cohort') +
  ylab('6WMD improvement (m)') +
  xlab('BMI group') +
  theme_bw() +
  stat_compare_means(method = "anova")



#for FC using differences in #makes sense, for walk difference not so much
crp_high_low_first_or_diagnostic_visit_v131Jan24 <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/crp_high_low_first_or_diagnostic_visit_v131Jan24.rds")

crp <- crp_high_low_first_or_diagnostic_visit_v131Jan24 %>% select(id, above_5)

teff <- merge(crp, treat2, by='id')
teff$above_5 <- as.factor(teff$above_5)
teff <- merge(teff, bmi, by='id')

mod1 <- lm(teff$walk_diff_meters ~ teff$bs_bmi)
mod2 <- lm(teff$fc_diff ~ teff$bs_bmi)

teff2 <- teff %>% filter(!is.na(bs_bmi) & !is.na(walk_diff_meters))
teff3 <- teff %>% filter(!is.na(bs_bmi) & !is.na(fc_diff))

teff2$resid_walk_diff <- mod1$residuals
teff3$resid_fc <- mod2$residuals

t.test(teff$fc_diff ~ teff$above_5)
t.test(teff$walk_diff_meters ~ teff$above_5)


hist(teff2$resid_walk_diff)
hist(teff3$resid_fc)

t.test(teff2$resid_walk_diff ~ teff2$above_5)
t.test(teff3$resid_fc ~ teff3$above_5)
#also NS for taking residuals for BMI

#regress out baseline walk
wl1 <- teff %>% filter(!is.na(ep_1_distance_meters_0) & !is.na(walk_diff_meters))
wl2 <- teff %>% filter(!is.na(ep_1_distance_meters_0) & !is.na(fc_diff))

m1 <- lm(wl1$walk_diff_meters ~ wl1$ep_1_distance_meters_0)
m2 <- lm(wl2$fc_diff ~ wl2$ep_1_distance_meters_0)

wl1$resid <- m1$residuals
wl2$resid <- m2$residuals

hist(wl1$resid)
hist(wl2$resid)

t.test(wl1$resid ~ wl1$above_5)
t.test(wl2$resid ~ wl2$above_5)

#test with removing BMI >40 
t.test(teff$ep_1_distance_meters_0, teff$ep_1_distance_meters_1, paired = TRUE)
t3 <- teff %>% filter(bs_bmi <=40)
t.test(t3$ep_1_distance_meters_0, t3$ep_1_distance_meters_1, paired = TRUE)

#also with BMI >40
t4 <- teff %>% filter(bs_bmi >40)
t.test(t4$ep_1_distance_meters_0, t4$ep_1_distance_meters_1, paired = TRUE)
#this difference is NS

#what if BMI >35
t4 <- teff %>% filter(bs_bmi >35)
t.test(t4$ep_1_distance_meters_0, t4$ep_1_distance_meters_1, paired = TRUE)
#this difference is sign
#===================================================================================================================
# Kawut paper: Gabler NM, French B, Strom BL, Palevsky HI, Taichman DB, Kawut SM, Halpern SD: 
#Validation of 6-minute walk distance as a surrogate end point in pulmonary arterial hypertension trials 
#Circulation 126 (3): 349-56,2012.

#identified threshold of 41.8m  for imrpovement
teff$sig_impove <- ifelse(teff$walk_diff_meters <41.8, 'no', 'yes')

kawut_teff <- teff %>% filter(!is.na(sig_impove))
#296 patients included in this
kawut_teff$sig_impove <- as.factor(kawut_teff$sig_impove)

#do some plotting
ggplot(kawut_teff, aes(x=sig_impove, y=bs_bmi)) +
  geom_boxplot()

#now get ROC curve for BMI for improvement
library(tidymodels)
model_fit <- multinom_reg() |> fit(sig_impove ~ bs_bmi, data=kawut_teff)

# show coefficients - shows no significance
tidy(model_fit, exponentiate = TRUE, conf.int = TRUE) |> 
  mutate_if(is.numeric, round, 4) |> 
  select(-std.error, -statistic)

glance(model_fit)
#deviance is not optimal....


pred_improve <- model_fit |> 
  augment(new_data = kawut_teff)

conf_mat(pred_improve, truth = sig_impove, estimate = .pred_class)
accuracy(pred_improve, truth = sig_impove, estimate = .pred_class)
#is accurate in about 50% of cases....... --> only predicts well for no improvement
roc_auc(pred_improve, truth = sig_impove, .pred_1)
roc_auc(pred_improve, truth = sig_impove, .pred_2) 


roc_1 <- roc(pred_improve, sig_impove, .pred_1)
coords(roc_1, "best", ret=c("threshold", "specificity", "sensitivity", '1-npv'))
# I dont think this sensitivity is correct as it makes no sense with th confusion matrix

plot(roc_1, print.auc=TRUE)
summary(roc_1)
