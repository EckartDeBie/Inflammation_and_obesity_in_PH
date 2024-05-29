#forest plot generation of linear models
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
library(grid)
library(forestploter)

#import the data

linear_models_aspire_adjusted <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/linear_models_aspire_adjusted.rds")
prepared_df_for_forest_models_cohort <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/prepared_df_for_forest_models_cohort.rds")

df1 <- linear_models_aspire_adjusted
df2 <- prepared_df_for_forest_models_cohort

df1$cohort <- 'ASPIRE'
df2$cohort <- 'UK cohort'

#rename dataframes to merge properly
df1$term <- ifelse(df1$term %in% c('log(bmi$level)', 'log(df2$level)'), 'CRP', 'BMI')
#also rename the tested variables
df1$variable_tested <- ifelse(df1$variable_tested == 'FVCp', 'FVC % predicted', df1$variable_tested)
df1$variable_tested <- ifelse(df1$variable_tested == 'BaseISWD', 'ISWD', df1$variable_tested)
df1$variable_tested <- ifelse(df1$variable_tested == 'CardiacOutput', 'Cardiac Output', df1$variable_tested)
df1$variable_tested <- ifelse(df1$variable_tested == 'CardiacIndex', 'Cardiac Index', df1$variable_tested)
df1$variable_tested <- ifelse(df1$variable_tested == 'bnp_level', 'BNP', df1$variable_tested)

df1$variable_tested <- paste(df1$variable_tested, df1$cohort, sep=" ")

#continue with df2
df2$term <- ifelse(df2$term == 'log(df$pref_crp)', 'CRP', 'BMI')
df2$variable_tested <- ifelse(df2$variable_tested == 'hb_pawp_m', 'PAWP', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'hb_pvr_calc', 'PVR', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'egfr_mdrd', 'eGFR', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'lf_fvc_pc', 'FVC % predicted', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'ep_1_distance_meters', '6MWD', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'pref_crp', 'CRP', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'hb_rap_m', 'mRAP', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'bs_bmi', 'BMI', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'hb_cardiac_index_value_1', 'Cardiac Index', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'hb_pap_m', 'mPAP', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'cbt_card_ntprobnp_ngpl', 'NT-proBNP', df2$variable_tested)
df2$variable_tested <- ifelse(df2$variable_tested == 'hb_cardiac_output_value_1', 'Cardiac Output', df2$variable_tested)

df2$variable_tested <- paste(df2$variable_tested, df2$cohort, sep=" ")

total_df <- rbind(df1, df2)
total_df$variable_tested <- as.factor(total_df$variable_tested)


total_df <- as.data.frame(total_df)

#now we can slowly plot
p <- ggplot(total_df, aes(y=total_df$variable_tested)) +
  geom_point(aes(x = total_df$estimate), shape=15, size=2) +
  geom_linerange(xmin = total_df$estimate-total_df$std.error, xmax = total_df$estimate+total_df$std.error) +
  geom_vline(xintercept = 0, linetype='dashed') +
  theme_classic() +
  labs(x='Standardised effect of variable', y='') +
  annotate("text", x = -.15, y = 16.5, label = "Negative association", size=3) +
  annotate("text", x = .15, y = 16.5, label = "Positive association", size=3) +
  ggtitle('Effect of BMI and CRP on outcome variables') +
  facet_wrap(~ term)

tdf <- total_df
tdf$` ` <- paste(rep(" ", 20), collapse = " ")


tdf <- separate(tdf, variable_tested, into = c("Outcome variable", "Cohort"), sep = " (?=[^ ]+$)")
tdf <- separate(tdf, `Outcome variable`, into = c("Outcome variable", "Cohort2"), sep = " (?=[^ ]+$)")

tdf$`Outcome variable` <- ifelse(tdf$`Outcome variable` == 'FVC %', 'FVC % predicted', tdf$`Outcome variable`)

tdf$ov <- tdf$`Outcome variable`

tdf <- tdf[order(tdf$term, tdf$ov),]

#make CSV to add in rows
write.csv2(tdf, file='standardised_effects_lm_for_plotting.csv')

#import the data
standardised_effects_lm_for_plotting <- read.csv2("~/PhD/Projects/CRP ~ survival and BMI/standardised_effects_lm_for_plotting.csv")
sdf <- standardised_effects_lm_for_plotting

#get high and low vars
sdf$high <- sdf$estimate + sdf$std.error
sdf$high <- round(sdf$high, digits = 4)
sdf$low <- sdf$estimate - sdf$std.error
sdf$low <- round(sdf$low, digits = 4)
sdf$estimate <- round(sdf$estimate, digits = 4)

#sdf$cohort <- ifelse(is.na(sdf$estimate), sdf$Outcome.variable, sdf$cohort)
#add in column
sdf$` ` <- paste(rep(" ", 20), collapse = " ")

#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
sdf$cohort <- ifelse(is.na(sdf$estimate), sdf$cohort, paste0("   ", sdf$cohort))

#make NAs blank
#sdf$p.value <- round(sdf$p.value, digits=6)
sdf$p.value <- ifelse(is.na(sdf$p.value), "", sdf$p.value)


#add CI column
sdf$`Estimate (95% CI)` <- ifelse(is.na(sdf$estimate), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   sdf$estimate, sdf$low, sdf$high))



tdfcrp <- sdf %>% filter(term == 'CRP')
tdfbmi <- sdf %>% filter(term == 'BMI')


p2 <- forest(tdfcrp[,c(10,22,23,6)],
            est = tdfcrp$estimate,
            lower = tdfcrp$low, 
            upper = tdfcrp$high,
            ci_column = 2,
            ref_line = 0,
            xlim = c(-0.4, 0.4),
            ticks_at = c(-0.4, -0.2, 0, 0.2, 0.4),
            xlab = 'Standardised effect on variable',
            title = 'Effect of CRP on clinical parameters',
            footnote = "Standardised effect of CRP on\nclinical parameters")
           

p3 <- forest(tdfbmi[,c(10,22,23,6)],
             est = tdfbmi$estimate,
             lower = tdfbmi$low, 
             upper = tdfbmi$high,
             ci_column = 2,
             ref_line = 0,
             xlab = 'Standardised effect on variable',
             xlim = c(-0.4, 0.4),
             ticks_at = c(-0.4, -0.2, 0, 0.2, 0.4),
             title='Effect of BMI on clinical parameters',
             footnote = "Standardised effect of BMI on\nclinical parameters")


# Print plot
plot(p2)
plot(p3)

library(gridExtra)
grid.arrange(p2, p3)

#=====================================================================
#now also get plot designed for paper
#import the data
standardised_effects_lm_for_plotting <- read.csv2("~/PhD/Projects/CRP ~ survival and BMI/standardised_effects_lm_for_plotting_reduced.csv")
sdf <- standardised_effects_lm_for_plotting

#get high and low vars
sdf$high <- sdf$estimate + sdf$std.error
sdf$high <- round(sdf$high, digits = 4)
sdf$low <- sdf$estimate - sdf$std.error
sdf$low <- round(sdf$low, digits = 4)
sdf$estimate <- round(sdf$estimate, digits = 4)

#sdf$cohort <- ifelse(is.na(sdf$estimate), sdf$Outcome.variable, sdf$cohort)
#add in column
sdf$` ` <- paste(rep(" ", 20), collapse = " ")

#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
sdf$cohort <- ifelse(is.na(sdf$estimate), sdf$cohort, paste0("   ", sdf$cohort))

#make NAs blank
#sdf$p.value <- round(sdf$p.value, digits=6)
sdf$p.value <- ifelse(is.na(sdf$p.value), "", sdf$p.value)


#add CI column
sdf$`Estimate (95% CI)` <- ifelse(is.na(sdf$estimate), "",
                                  sprintf("%.2f (%.2f to %.2f)",
                                          sdf$estimate, sdf$low, sdf$high))



tdfcrp <- sdf %>% filter(term == 'CRP')
tdfbmi <- sdf %>% filter(term == 'BMI')


p2 <- forest(tdfcrp[,c(10,15)],
             est = tdfcrp$estimate,
             lower = tdfcrp$low, 
             upper = tdfcrp$high,
             ci_column = 2,
             ref_line = 0,
             xlim = c(-0.4, 0.4),
             ticks_at = c(-0.4, -0.2, 0, 0.2, 0.4),
             xlab = 'Standardised effect on variable',
             title = 'Effect of CRP on clinical parameters')
         
p3 <- forest(tdfbmi[,c(10,15)],
             est = tdfbmi$estimate,
             lower = tdfbmi$low, 
             upper = tdfbmi$high,
             ci_column = 2,
             ref_line = 0,
             xlab = 'Standardised effect on variable',
             xlim = c(-0.4, 0.4),
             ticks_at = c(-0.4, -0.2, 0, 0.2, 0.4),
             title='Effect of BMI on clinical parameters')
             



#====================================================================
#now do the forest plots for the coxPH

#import the data
coefficient_cox_ph_CRP_groups_aspire_v1_7May <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/coefficient_cox_ph_CRP_groups_aspire_v1_7May.rds")
survival_10_years_post_diagnosis_coxPH_crp_cohort_v1_7May <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/survival_10_years_post_diagnosis_coxPH_crp_cohort_v1_7May.rds")

#save files in CSV and prepare for forest plot
write.csv2(coefficient_cox_ph_CRP_groups_aspire_v1_7May, 'aspire_cox_crp.csv')
write.csv2(survival_10_years_post_diagnosis_coxPH_crp_cohort_v1_7May,'cohort_cox_crp.csv')

#import the formatted data
Cohort_ASPIRE_CRP_coxPH <- read_excel("CRP ~ survival and BMI/Cohort_ASPIRE_CRP_coxPH.xlsx")
dfpl <- Cohort_ASPIRE_CRP_coxPH


dfpl$` ` <- paste(rep(" ", 20), collapse = " ")

#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
dfpl$Variable <- ifelse(is.na(dfpl$lower), dfpl$Variable, paste0("   ", dfpl$Variable))

#make NAs blank
#sdf$p.value <- round(sdf$p.value, digits=6)
dfpl$`p-value` <- ifelse(is.na(dfpl$`p-value`), "", dfpl$`p-value`)


#add CI column
dfpl$`HR (95% CI)` <- ifelse(is.na(dfpl$HR), "",
                                  sprintf("%.2f (%.2f to %.2f)",
                                          dfpl$HR, dfpl$lower, dfpl$upper))

p3 <- forest(dfpl[,c(1,6,7,5)],
             est = dfpl$HR,
             lower = dfpl$lower, 
             upper = dfpl$upper,
             ci_column = 2,
             ref_line = 1,
             xlab = 'Hazard ratio with 95% CI',
             xlim = c(0.25, 2),
             ticks_at = c(0.25, 0.5, 1, 1.5, 2),
             title='Cox-proportional Hazard model for IPAH/HPAH and ASPIRE cohorts',
             footnote = "HR wiht 95% CI")


#do per PH group
#load in data
coefficient_cox_ph_CRP_groups_aspire_v1_7May <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/coefficient_cox_ph_CRP_groups_aspire_v1_7May.rds")
summary_PH2_cox_aspire <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/summary_PH2_cox_aspire.rds")
summary_PH3_cox_aspire <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/summary_PH3_cox_aspire.rds")
coefPAH <- as.data.frame(summary_PAH_cox_aspire$coefficients)
summary_PH4_cox_aspire <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/summary_PH4_cox_aspire.rds")
summary_PH5_cox_aspire <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/CoxPH coefficients/summary_PH5_cox_aspire.rds")
ph2 <- summary_PH2_cox_aspire
ph3 <- summary_PH3_cox_aspire
ph4 <- summary_PH4_cox_aspire
ph5 <- summary_PH5_cox_aspire

coefPAH$group <- 'Group 1 PAH'
ph2$group <- 'Group 2 PH'
ph3$group <- 'Group 3 PH'
ph4$group <- 'Group 4 PH'
ph5$group <- 'Group 5 PH'

ph <- rbind(coefPAH, ph2, ph3,ph4, ph5)
ph$upper_ci <- ph$`exp(coef)` + (1.96*ph$`se(coef)`)
ph$lower_ci <- ph$`exp(coef)` - (1.96*ph$`se(coef)`)
ph$variable <- rownames(ph)

ph <- ph[order(ph$variable),]
#make csv
write.csv2(ph, file='coefficients_PH_subgroups_Cox.csv')

#import edited csv file
coefficients_PH_subgroups_Cox <- read.csv("~/PhD/Projects/CRP ~ survival and BMI/coefficients_PH_subgroups_Cox.csv", sep=";")

df2 <- coefficients_PH_subgroups_Cox

df2$` ` <- paste(rep(" ", 20), collapse = " ")

#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
df2$group <- ifelse(is.na(df2$HR), df2$group, paste0("   ", df2$group))

#make NAs blank
df2$P.value <- ifelse(is.na(df2$P.value), "", df2$P.value)


#add CI column
df2$`HR (95% CI)` <- ifelse(is.na(df2$HR), "",
                             sprintf("%.2f (%.2f to %.2f)",
                                     df2$HR, df2$lower_ci, df2$upper_ci))

p3 <- forest(df2[,c(4,8,9,3)],
             est = df2$HR,
             lower = df2$lower_ci, 
             upper = df2$upper_ci,
             ci_column = 2,
             ref_line = 1,
             xlab = 'Hazard ratio with 95% CI',
             xlim = c(0.7, 2.5),
             ticks_at = c(0.7, 1, 1.5, 2, 2.5),
             title='Cox-proportional Hazard model for PH subgroup in ASPIRE',
             footnote = "HR wiht 95% CI")

#==============================================================================
#EdB continued here on 15-05-2024
#scaled forest plots

#import the data
coefficients_SCALED_cox_per_CRP_cohort <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/SCALED coefficients coxPH/coefficients_SCALED_cox_per_CRP_cohort.rds")
SCALED_coefficients_coxPHallPH_per_CRP <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/SCALED coefficients coxPH/SCALED_coefficients_coxPHallPH_per_CRP.rds")

c1 <- coefficients_SCALED_cox_per_CRP_cohort
c1$group <- 'UK cohort'
c2 <- SCALED_coefficients_coxPHallPH_per_CRP
c2$group <- 'ASPIRE'


c1$lower <- c1$`exp(coef)` - (c1$`se(coef)` * 1.96)
c1$upper <- c1$`exp(coef)` + (c1$`se(coef)` * 1.96)

c2$lower <- c2$`exp(coef)` - (c2$`se(coef)` * 1.96)
c2$upper <- c2$`exp(coef)` + (c2$`se(coef)` * 1.96)
  
write.csv2(c1, file='cohort_scaled_coef_crp.csv')
write.csv2(c2, file='aspire_scaled_coef_crp.csv')


#now plot!
SCALED_crp_allPH_coefs <- read.csv2("~/PhD/Projects/CRP ~ survival and BMI/SCALED_crp_allPH_coefs.csv")
df2 <- SCALED_crp_allPH_coefs

names(df2) <- c('x', 'no', 'HR', 'se', 'z', 'P.value', 'group', 'lower_ci', 'upper_ci')

df2$` ` <- paste(rep(" ", 20), collapse = " ")


#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
df2$group <- ifelse(is.na(df2$HR), df2$group, paste0("     ", df2$group))

#add CI column
df2$`HR (95% CI)` <- ifelse(is.na(df2$HR), "",
                            sprintf("%.2f (%.2f to %.2f)",
                                    df2$HR, df2$lower_ci, df2$upper_ci))

#round pvals and otherwise give them another number
df2$P.value <- round(df2$P.value, digits=5)

#make NAs blank
df2$P.value <- ifelse(df2$P.value < 0.00001, '<0.00001', df2$P.value)
df2$P.value <- ifelse(df2$P.value == 1e-05, '<0.00001', df2$P.value)
df2$P.value <- ifelse(df2$P.value == 2e-05, '<0.00001', df2$P.value)

df2$`P-value` <- df2$P.value

#make P-value blank if NA
df2$`P-value` <- ifelse(is.na(df2$`P-value`), "", df2$`P-value`)

p3 <- forest(df2[,c(7,10,11,12)],
             est = df2$HR,
             lower = df2$lower_ci, 
             upper = df2$upper_ci,
             ci_column = 2,
             ref_line = 1,
             xlab = 'Hazard ratio with 95% CI',
             xlim = c(0, 2),
             ticks_at = c(0, 0.5, 1, 1.5, 2),
             title='Cox-proportional Hazard model for the effect of diagnostic CRP in PH',
             footnote = "HR wiht 95% CI")


#now repeat this for the CRP effects per PH subtype in ASPIRE
coefficients_SCALED_effects_cox_crp_per_PH_subgroup <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/SCALED coefficients coxPH/coefficients_SCALED_effects_cox_crp_per_PH_subgroup.rds")

coef <- coefficients_SCALED_effects_cox_crp_per_PH_subgroup
coef$lower_ci <- coef$`exp(coef)` - (coef$`se(coef)`*1.96)
coef$upper_ci <- coef$`exp(coef)` + (coef$`se(coef)`*1.96)
coef$var <- rownames(coef)
coef <- coef[order(coef$var),]

write.csv2(coef, file='coef_SCALED_per_PH_group_crp.csv')

#now import the data
coef_SCALED_per_PH_group_crp <- read.csv2("~/PhD/Projects/CRP ~ survival and BMI/coef_SCALED_per_PH_group_crp.csv")
df2 <- coef_SCALED_per_PH_group_crp


names(df2) <- c('x', 'no', 'HR', 'se', 'z', 'P.value', 'group', 'lower_ci', 'upper_ci', 'var')

df2$` ` <- paste(rep(" ", 20), collapse = " ")


#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
df2$group <- ifelse(is.na(df2$HR), df2$group, paste0("     ", df2$group))

#add CI column
df2$`HR (95% CI)` <- ifelse(is.na(df2$HR), "",
                            sprintf("%.2f (%.2f to %.2f)",
                                    df2$HR, df2$lower_ci, df2$upper_ci))

#round pvals and otherwise give them another number
df2$P.value <- round(df2$P.value, digits=5)

#make NAs blank
df2$P.value <- ifelse(df2$P.value < 0.0001, '<0.0001', df2$P.value)
df2$P.value <- ifelse(df2$P.value == 1e-05, '<0.0001', df2$P.value)
df2$P.value <- ifelse(df2$P.value == 2e-05, '<0.0001', df2$P.value)
df2$P.value <- ifelse(df2$P.value %in% c(9e-05, 7e-05, 2e-04), '<0.0001', df2$P.value)



df2$`P-value` <- df2$P.value

#make P-value blank if NA
df2$`P-value` <- ifelse(is.na(df2$`P-value`), "", df2$`P-value`)

p3 <- forest(df2[,c(7,11,12,13)],
             est = df2$HR,
             lower = df2$lower_ci, 
             upper = df2$upper_ci,
             ci_column = 2,
             ref_line = 1,
             xlab = 'Hazard ratio with 95% CI',
             xlim = c(0, 2),
             ticks_at = c(0, 0.5, 1, 1.5, 2),
             title='Cox-proportional Hazard model per PH subgroup',
             footnote = "HR wiht 95% CI")

#===========================================================================
#BMI CoxPH
#get the CoxPH data

SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_v1_3May24 <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/SCALED coefficients coxPH/SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_v1_3May24.rds")

sccoef <- SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_v1_3May24
sccoef$lower <- sccoef$`exp(coef)` - (sccoef$`se(coef)`*1.96)
sccoef$upper <- sccoef$`exp(coef)` + (sccoef$`se(coef)`*1.96)

#make CSV and adapt by hand
write.csv2(sccoef, file='coefs_BMI_for_plot.csv')


#import the data
coefs_BMI_for_plot <- read.csv2("~/PhD/Projects/CRP ~ survival and BMI/coefs_BMI_for_plot.csv")

c1 <- coefs_BMI_for_plot
df2 <- c1


df2$` ` <- paste(rep(" ", 20), collapse = " ")


#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
df2$Group <- ifelse(is.na(df2$HR), df2$Group, paste0("     ", df2$Group))

df2$lower <- round(df2$lower, digits=4)
df2$upper <- round(df2$upper, digits=4)

#add CI column
df2$`HR (95% CI)` <- ifelse(is.na(df2$HR), "",
                            sprintf("%.2f (%.2f to %.2f)",
                                    df2$HR, df2$lower, df2$upper))

#round pvals and otherwise give them another number
df2$P.value <- round(df2$P.val, digits=5)

#make NAs blank
df2$P.value <- ifelse(df2$P.value < 0.0001, '<0.0001', df2$P.value)
df2$P.value <- ifelse(df2$P.value %in% c(9e-04, 7e-05, 2e-04), '<0.0001', df2$P.value)



df2$`P-value` <- df2$P.value

#make P-value blank if NA
df2$`P-value` <- ifelse(is.na(df2$`P-value`), "", df2$`P-value`)

p3 <- forest(df2[,c(6,11,12,14)],
             est = df2$HR,
             lower = df2$lower, 
             upper = df2$upper,
             ci_column = 2,
             ref_line = 1,
             xlab = 'Hazard ratio with 95% CI',
             xlim = c(0, 6),
             ticks_at = c(0, 1, 2, 3, 4, 5,6),
             title='Cox-proportional Hazard model for BMI group and covariates',
             footnote = "HR wiht 95% CI")
#================================================
#do cofficients with the 1 year removed plot
#import the data
SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_1yr_removed_v1_7May24 <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/SCALED coefficients coxPH/SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_1yr_removed_v1_7May24.rds")
coef1 <- SCALED_Coefficients_ASPIRE_Cox_PH_per_BMI_group_1yr_removed_v1_7May24
coef1$upper <- coef1$`exp(coef)` + (coef1$`se(coef)` * 1.96)
coef1$lower <- coef1$`exp(coef)` - (coef1$`se(coef)` * 1.96)

write.csv2(coef1, file='coefs_to_plot_BMI_1yr_removed.csv')
#manually changed layout --> now import
coefs_to_plot_BMI_1yr_removed <- read.csv2("~/PhD/Projects/CRP ~ survival and BMI/coefs_to_plot_BMI_1yr_removed.csv")

df2 <- coefs_to_plot_BMI_1yr_removed


df2$` ` <- paste(rep(" ", 20), collapse = " ")


#indent cohort (https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)
df2$Group <- ifelse(is.na(df2$HR), df2$Group, paste0("     ", df2$Group))

df2$lower <- round(df2$lower, digits=4)
df2$upper <- round(df2$upper, digits=4)

#add CI column
df2$`HR (95% CI)` <- ifelse(is.na(df2$HR), "",
                            sprintf("%.2f (%.2f to %.2f)",
                                    df2$HR, df2$lower, df2$upper))

#round pvals and otherwise give them another number
df2$P.value <- round(df2$P.val, digits=5)

#make NAs blank
df2$P.value <- ifelse(df2$P.value < 0.0001, '<0.0001', df2$P.value)

df2$`P-value` <- df2$P.value

#make P-value blank if NA
df2$`P-value` <- ifelse(is.na(df2$`P-value`), "", df2$`P-value`)

p3 <- forest(df2[,c(6,7,8,10)],
             est = df2$HR,
             lower = df2$lower, 
             upper = df2$upper,
             ci_column = 2,
             ref_line = 1,
             xlab = 'Hazard ratio with 95% CI',
             xlim = c(0, 2),
             ticks_at = c(0, 0.5, 1, 1.5, 2),
             title='Cox-proportional Hazard model for BMI; 1st year of survival removed',
             footnote = "HR wiht 95% CI")
