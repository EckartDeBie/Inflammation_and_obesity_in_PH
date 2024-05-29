#odds ratio's for smoking and comorbid diseases in BMI groups
#CRP modelling without clustering
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
getwd()

#load packages
library("tidyverse")

#load the data
v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")
cleaned_comorbidity_dataV2_20Dec2023 <- readRDS("C:/Users/location/cleaned_comorbidity_dataV2_20Dec2023.rds")
query_smoking_dump_10_7_23 <- read.delim("~/PhD/Projects/CRP ~ survival and BMI/query_smoking_dump_10_7_23.txt")

df1 <- v4_clean_clinical_data_first_visit_15Dec23 %>% select(id, bs_bmi)
#remove underweight patients
df1 <- df1 %>% filter(bs_bmi >=18.5)
df1$weight <- ifelse(df1$bs_bmi >=30, 'obese', 'normal weight or overweight')
df1$weight <- as.factor(df1$weight)

com <- cleaned_comorbidity_dataV2_20Dec2023
com$htn <-  ifelse(com$cfh_comorbid_disease_diagnosis == 'hypertension', 'yes', 'no')
com$dm2 <- ifelse(com$cfh_comorbid_disease_diagnosis == 'DM2', 'yes', 'no')
com$hypothyr <- ifelse(com$cfh_comorbid_disease_diagnosis == 'hypothyroidism', 'yes', 'no')
com$dyslip <- ifelse(com$cfh_comorbid_disease_diagnosis == 'dyslipidaemia', 'yes', 'no')
com$iheart <- ifelse(com$cfh_comorbid_disease_diagnosis == 'ischaemic_heart_disease', 'yes', 'no')
com$copd <- ifelse(com$cfh_comorbid_disease_diagnosis == 'copd', 'yes', 'no')
com$asthma <- ifelse(com$cfh_comorbid_disease_diagnosis == 'asthma', 'yes', 'no')
com$pe <- ifelse(com$cfh_comorbid_disease_diagnosis == 'pulmonary embolism', 'yes', 'no')
com$ckd <- ifelse(com$cfh_comorbid_disease_diagnosis == 'CKD', 'yes', 'no')
com$obese <- ifelse(com$cfh_comorbid_disease_diagnosis == 'obesity', 'yes', 'no')
com$mi <- ifelse(com$cfh_comorbid_disease_diagnosis == 'myocardial infarction', 'yes', 'no')
com$osa <- ifelse(com$cfh_comorbid_disease_diagnosis == 'osa', 'yes', 'no')
com$iron <- ifelse(com$cfh_comorbid_disease_diagnosis == 'iron_deficiency', 'yes', 'no')



smoke2 <- query_smoking_dump_10_7_23 %>% select(study_subject_oid, cfh_smoking_history, cfh_smoking_current)
smoke2 <- unique(smoke2)
smoke2$cfh_smoking_history <- ifelse(smoke2$cfh_smoking_history == 'NULL', 'no', smoke2$cfh_smoking_history)
smoke2$cfh_smoking_history <- ifelse(smoke2$cfh_smoking_history == 'UNK', NA, smoke2$cfh_smoking_history)
smoke2$cfh_smoking_history <- as.factor(smoke2$cfh_smoking_history)

smoke2$smoke_now <- ifelse(smoke2$cfh_smoking_history == 'yes', 'past_smoker', 'non_smoker')
smoke2$smoke_now <- ifelse(smoke2$cfh_smoking_current == 'yes', 'current_smoker', smoke2$smoke_now)
smoke2$smoke_now <- as.factor(smoke2$smoke_now)

smokers <- smoke2 %>% select(study_subject_oid, smoke_now)

com2 <- com[,c(1,36:45,47:49)]
com3 <- merge(com2, smokers, by.x='id', by.y='study_subject_oid')

#now merge with BMI data
com4 <- merge(com3, df1)
com4[,2:14] <- lapply(com4[,2:14], as.factor)

com5 <- com4 %>% group_by(id) %>% summarise(across(everything(), ~ if("yes" %in% .x) "yes" else "no"))
#this has worked
com5 <- unique(com5)
com5 <- com5[,-c(15:17)]
#select other identifiers from com4

com6 <- merge(com5, com4[,c(1,15:17)])
com6 <- unique(com6)
#1001 patients retained

#now get tables
library(Publish)
a1 <- univariateTable(weight ~ htn + dm2 + hypothyr + dyslip + iheart + copd + asthma + pe + ckd + mi + osa + iron + smoke_now, data=com6, show.totals = FALSE, column.percent = FALSE, compare.groups = FALSE)
a2 <- summary(a1)

df_cleaned <- a2 %>% mutate(across(where(is.character), ~gsub("\\s*\\([^\\)]+\\)", "", .x)))
df_cleaned$`normal weight or overweight (n=613)` <- as.numeric(df_cleaned$`normal weight or overweight (n=613)`)
df_cleaned$`obese (n=388)` <- as.numeric(df_cleaned$`obese (n=388)`)

#calculate OR in excel
write.csv2(df_cleaned, 'number_of_comorbidities_in_obesity_groups.csv')

#import the calculated ORs with SE to plot
Obesity_odds_ratios <- read_excel("CRP ~ survival and BMI/Obesity_odds_ratios.xlsx")

library(grid)
library(forestploter)

OR <- Obesity_odds_ratios
OR$OR <- round(OR$OR, digits = 4)
OR$upper <- round(OR$upper, digits =4)
OR$lower <- round(OR$lower, digits=4)

names(OR) <- c('Comorbid condition', 'OR', '95% CI upper bound', '95% CI lower bound')

OR$` ` <- paste(rep(" ", 20), collapse = "        ")


p <- forest(OR[,c(1,5,2,3,4)],
            est = OR$OR,
            lower = OR$`95% CI lower bound`, 
            upper = OR$`95% CI upper bound`,
            ci_column = 2,
            ref_line = 1,
            xlab = 'Odds ratio with 95% CI',
            xlim = c(0, 12),
            title = 'Odds ratios for comorbid diseases for obesity at diagnosis based on BMI',
            ticks_at = c(0, 2, 4, 6, 8, 10, 12),
            footnote = "Odds ratio for disease based on\ndiagnostic BMI")
