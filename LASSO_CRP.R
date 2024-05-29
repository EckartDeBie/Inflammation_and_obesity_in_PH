#LASSO model for REVEAL risk score
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
getwd()

#load packages
library("lubridate")
library("ggsci")
library("magrittr")
library("tidyverse")
#load in the data
v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")
v4_clean_clinical_data_first_visit_15Dec23$SVI <- v4_clean_clinical_data_first_visit_15Dec23$hb_cardiac_output_value_1/v4_clean_clinical_data_first_visit_15Dec23$hb_heart_rate
v4_clean_clinical_data_first_visit_15Dec23$dlco <- ((v4_clean_clinical_data_first_visit_15Dec23$lf_kco_pc/100)*(v4_clean_clinical_data_first_visit_15Dec23$lf_va_pc/100)*100)

load("C:/Users/Gebruiker/Downloads/data_clean (3).RData")
centre <- data.clean.cohort %>% select('id', 'centre')
centre <- centre %>% filter(centre %in% c('Glasgow', 'Sheffield', 'Great Ormond Street', 'Lincoln', 'Papworth', 'Royal Brompton', 'Royal United Hospital Bath', 'Imperial and Hammersmith', 'Newcastle Freeman', 'Royal Free'))
centre <- unique(centre)

df <- merge(v4_clean_clinical_data_first_visit_15Dec23, centre, by='id')
load("C:/Users/Gebruiker/Downloads/data_checked (4).RData")
sel <- data.checked %>% select(id, cbt_card_bnp_ngpl)
da <- merge(df, sel, by='id')

df <- da %>% select(id, DOB, age_diagnosis, sub_date, sub_cause, sex, functional_class, ep_1_distance_meters, ep_1_type_6mwt, hb_pap_s, hb_rap_m, hb_cardiac_index_value_1, hb_sv_o2, SVI, ec_right_atrial_area, ec_percicardial_effusion, ec_tricuspid_apse, cbt_card_ntprobnp_ngpl, cbt_card_bnp_ngpl)

mort_df <- df
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
mort_df <- mort_df %>% filter(!is.na(surv_time))
#now censor the survival times
mort_df$event <- ifelse(mort_df$surv_time >=10, 0, mort_df$event)
mort_df$surv_time <- ifelse(mort_df$surv_time >=10, 10, mort_df$surv_time)
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)

surv <- mort_df %>% select(id, diagnosis_date)

comlete_df <- merge(surv, df, by='id')

complete_df <- comlete_df[,-c(3)]



#get REVEAL scores per col
add_phenotype_reveal_risk_group <- function(dt) {
  n_patients <- nrow(dt)
  reveal <- matrix(NA, nrow = n_patients, ncol = 13)
  
  # Loop to calculate scores based on conditions
  for (i in seq_len(n_patients)) {
    reveal[i, 1] = ifelse(dt$diagnosis_verified[i] == "IPAH", 0, 2)
    reveal[i, 2] = ifelse(is.na(dt$age_diagnosis[i]) & dt$sex[i] == "male", NA, 
                          ifelse(dt$sex[i] == "male" & dt$age_diagnosis[i] <= 60 | dt$sex[i] == "female", 0, 2))
    reveal[i,3]=ifelse(is.na(dt$egfr_mdrd)[i], NA, ifelse(dt$egfr_mdrd[i] < 60, 1, 0))
    reveal[i,4]=ifelse(is.na(dt$functional_class)[i], NA, ifelse(dt$functional_class[i] == "2", 0, ifelse(dt$functional_class[i] == "1", -2, ifelse(dt$functional_class[i] == "3", 1, 2))))
    reveal[i,5]=ifelse(is.na(dt$cfe_bp_systolic)[i], NA, ifelse(dt$cfe_bp_systolic[i] < 110, 1, 0))
    reveal[i,6]=ifelse(is.na(dt$cfe_heart_rate)[i], NA, ifelse(dt$cfe_heart_rate[i] > 96, 1, 0))
    reveal[i,7]=ifelse(is.na(dt$ep_1_distance_meters)[i], NA, ifelse(dt$ep_1_distance_meters[i] < 165, 1, ifelse(dt$ep_1_distance_meters[i] >= 440, -1, 0)))
    reveal[i,8]=ifelse(is.na(dt$ec_percicardial_effusion)[i], NA, ifelse(dt$ec_percicardial_effusion[i] == "absent", 0, 1))
    reveal[i,9]=ifelse(is.na(dt$dlco)[i], NA, ifelse(dt$dlco[i] < 40, 1, 0))
    reveal[i,10]=ifelse(is.na(dt$hb_rap_m)[i], NA, ifelse(dt$hb_rap_m[i] > 20, 1, 0))
    reveal[i,11]=ifelse(is.na(dt$hb_pvr_calc)[i], NA, ifelse(dt$pvr[i] <400, -1, 0))
    reveal[i,12]=ifelse(is.na(dt$cbt_card_ntprobnp_ngpl)[i], NA, ifelse(dt$cbt_card_ntprobnp_ngpl[i] <300, -2, ifelse(dt$cbt_card_ntprobnp_ngpl[i] >1100, 2, 0)))
    reveal[i,13]=ifelse(is.na(dt$cbt_card_bnp_ngpl)[i], NA, ifelse(dt$cbt_card_bnp_ngpl[i] < 50 | dt$cbt_card_bnp_ngpl[i] <200, -2, ifelse(dt$cbt_card_bnp_ngpl[i] > 200 | dt$cbt_card_bnp_ngpl[i] >=800, 2, 0)))
  }
  
  # Convert matrix to data frame and add column names
  reveal_score <- as.data.frame(reveal)
  colnames(reveal_score) <- c("diagnosis", "sex", "renal_failure", "WHO", "SBP", "HR", "walk",
                              "pericard_effusion", "DLCO", "mRAP", "PVR", "NTproBNP", "BNP")
  
  # Calculate total NAs and total score
  reveal_score$total_na <- rowSums(is.na(reveal_score))
  reveal_score$total_score <- rowSums(reveal_score, na.rm = TRUE, dims = 1) + 6
  
  # Add patient IDs
  reveal_score$id <- dt$id
  
  # Return the merged dataset
  dt <- merge(dt, reveal_score, by = "id", sort = FALSE)
  return(dt)
}

da2 <- add_phenotype_reveal_risk_group(da)
#use same col for BNP/NTproBNP
da2$NTproBNP <- ifelse(is.na(da2$NTproBNP), da2$BNP, da2$NTproBNP)

#add in score for CRP and for BMI
da2$reveal_plus_CRP <- NA
da2$hl_crp <-  ifelse(da2$cbt_inflammation_crp_mgpl > 5, 'high', 'low')
da2$reveal_plus_CRP <-  ifelse(da2$hl_crp == 'high', 2, -2)

da2$reveal_plus_BMI <- NA
da2$reveal_plus_BMI <- ifelse(da2$bs_bmi <18.5, + 2, da2$reveal_plus_BMI)
da2$reveal_plus_BMI <- ifelse(da2$bs_bmi >=25 & da2$bs_bmi <30, -1, da2$reveal_plus_BMI)
da2$reveal_plus_BMI <- ifelse(da2$bs_bmi >=30, -2, da2$reveal_plus_BMI)

#select the relevant cols 
df_forlas <- da2[,c(1,10, 318:329, 333,335, 331)]
#only select UK patients for survival
df_forlas <- df_forlas %>% filter(id %in% centre$id)

#remove  >= 5 points missing from score
df_forlas <- df_forlas %>% filter(total_na <5)
#961 left
#also only want adult patients
df_forlas <- df_forlas %>% filter(age_diagnosis >=18)
#895 patients left for an analysis

#merge with mortality data
df2 <- mort_df[,c(1,22,24)]

df_forlas <- merge(df_forlas[,c(1,3:16)], df2, by='id')
#888 kept in

library(glmnet)
#==============================================================
#now do LASSO to predict survival time (censored at 10 years)
#now select a subset of the df to work with
df3 <- na.omit(df_forlas)
#retains a too small number to just work with....
#134 patients are kept in --> need to use as cannot have NAs

rownames(df3) <- df3$id
#select 75% randomly and keep the other 25% for validation
set.seed(123)
sample_df <- df3[sample(nrow(df3), 100),]
#now we've got a random subset

#get response variable 
resp <- sample_df$surv_time
#get predictor variables
pred <- data.matrix(sample_df[,c(2:15)])
#==========================================================
#try to explain life status
#=========================================================
resp <- as.character(sample_df$event)
resp <- as.factor(sample_df$event)

#run the lasso model
lasso <- cv.glmnet(pred, resp, alpha = 1, nfolds = 10, family='multinomial')

#use K-fold cross validation for lambda
best_lambda <- lasso$lambda.min

plot(lasso)

#now run glmnet
mod <- glmnet(pred, resp, alpha = 1, lambda = best_lambda, family = 'multinomial')
coef(mod)

#check how often this works
val2 <- df_forlas %>% filter(!id %in% sample_df$id)
rownames(val2) <- val2$id
val <- data.matrix(na.omit(val2[,c(2:15)]))

predval <- predict(mod, s=best_lambda, newx = val, type = 'class')
pv <- as.data.frame(predval)
pv$id <- rownames(val)

pv <- merge(pv, val2[,c(1,17)])
names(pv)[[2]] <- 'lasso_pred'
pv$lasso_pred <- as.factor(pv$lasso_pred)
#see if the labels correspond
pv$same <- ifelse(pv$lasso_pred == pv$event, 'same', 'not same')
summary(as.factor(pv$same))
#24/34 are accurate; 70%...... This is with all the data

#now see what happens to ROC curves with less information --> only taking sex, egfr, sbp, walk distance, pericardial effusion, DLCO, NTproBNP and CRP

df_rev <- df_forlas
#replace all NA with 0 (as already adults + <5 missing for REVEAL score)
df_rev[is.na(df_rev)] <- 0
#now add
df_rev$REVEAL2 <- df_rev$diagnosis + df_rev$sex.y + df_rev$renal_failure + df_rev$WHO + df_rev$SBP + df_rev$SBP + df_rev$HR + df_rev$walk + df_rev$pericard_effusion + df_rev$DLCO + df_rev$mRAP + df_rev$PVR + df_rev$NTproBNP
df_rev$reveal_crp <- df_rev$REVEAL2 + df_rev$reveal_plus_CRP
df_rev$reveal_bmi <- df_rev$REVEAL2 + df_rev$reveal_plus_BMI
df_rev$reveal_crp_bmi <- df_rev$REVEAL2 + df_rev$reveal_plus_CRP +df_rev$reveal_plus_BMI
df_rev$Noninvasive_reveal <- df_rev$sex.y + df_rev$renal_failure + df_rev$SBP + df_rev$walk + df_rev$pericard_effusion + df_rev$DLCO + df_rev$NTproBNP 
df_rev$Noninvasive_reveal_crp <- df_rev$sex.y + df_rev$renal_failure + df_rev$SBP + df_rev$walk + df_rev$pericard_effusion + df_rev$DLCO + df_rev$NTproBNP + df_rev$reveal_plus_CRP
df_rev$Noninvasive_reveal_crp_bmi <- df_rev$sex.y + df_rev$renal_failure + df_rev$SBP + df_rev$walk + df_rev$pericard_effusion + df_rev$DLCO + df_rev$NTproBNP + df_rev$reveal_plus_CRP + df_rev$reveal_plus_BMI


library(survivalROC)
library(timeROC)

roc1 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$REVEAL2, predict.time = 10, method = "KM")
roc2 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$reveal_crp, predict.time = 10, method = "KM")
roc3 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$reveal_bmi, predict.time = 10, method = "KM")
roc4 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$reveal_crp_bmi, predict.time = 10, method = "KM")
roc5 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$Noninvasive_reveal, predict.time = 10, method = "KM")
roc6 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$Noninvasive_reveal_crp, predict.time = 10, method = "KM")
roc7 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$Noninvasive_reveal_crp_bmi, predict.time = 10, method = "KM")


#now plot
#now plot
plot(roc1$FP, roc1$TP, type = "l", col = "yellow", xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC Curves Comparison of REVEAL score 10 years post diagnosis UK cohort")
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


#as reveal normally only calculates 12 month survival
roc1 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$REVEAL2, predict.time = 1, method = "KM")
roc2 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$reveal_crp, predict.time = 1, method = "KM")
roc3 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$reveal_bmi, predict.time = 1, method = "KM")
roc4 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$reveal_crp_bmi, predict.time = 1, method = "KM")
roc5 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$Noninvasive_reveal, predict.time = 1, method = "KM")
roc6 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$Noninvasive_reveal_crp, predict.time = 1, method = "KM")
roc7 <- survivalROC(Stime = df_rev$surv_time, status = df_rev$event, marker = df_rev$Noninvasive_reveal_crp_bmi, predict.time = 1, method = "KM")


#now plot
#now plot
plot(roc1$FP, roc1$TP, type = "l", col = "red", xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC Curves Comparison of REVEAL score 1 year post diagnosis UK cohort")
#lines(roc2$FP, roc2$TP, col = "black")
#lines(roc3$FP, roc3$TP, col = "yellow")
#lines(roc4$FP, roc4$TP, col = "green")
lines(roc5$FP, roc5$TP, col = "blue")
#lines(roc6$FP, roc6$TP, col = "cyan")
lines(roc7$FP, roc7$TP, col = "darkgreen")

legend("bottomright", legend = c(paste("REVEAL 2.0, AUC:", round(roc1$AUC, 2)),
                                 paste("Noninvasive REVEAL, AUC:", round(roc5$AUC, 2)),
                                 paste("Noninvasive REVEAL + CRP & BMI, AUC:", round(roc7$AUC, 2))),
       col = c("red", "blue", "darkgreen"), lty = 1, cex=0.8)


