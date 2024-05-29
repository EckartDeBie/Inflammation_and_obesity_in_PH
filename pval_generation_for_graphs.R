#pval generation of clinical differences CRP clusters
#CRP modelling without clustering
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
getwd()

#load packages
library("tidyverse")

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")
log_transformed_crp_levels_clusters <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/log_transformed_crp_levels_clusters.rds")


clust <- log_transformed_crp_levels_clusters
names(clust) <- c('id', 'cluster')
clust$cluster <- as.factor(clust$cluster)

df <- merge(v4_clean_clinical_data_first_visit_15Dec23, clust, by='id')

wilcox.test(df$bs_bmi ~ df$cluster)
wilcox.test(df$hb_pawp_m ~ df$cluster)
wilcox.test(df$hb_rap_m ~ df$cluster)
wilcox.test(df$egfr_mdrd ~ df$cluster)
t.test(df$lf_fvc_pc ~ df$cluster)

df2 <- df %>% filter(ep_1_type_6mwt == 'corridor')
wilcox.test(df2$ep_1_distance_meters ~ df2$cluster)


#==============================================================
#repeat analysis for high vs low groups with initial CRP
crp_high_low_first_or_diagnostic_visit_v131Jan24 <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/crp_high_low_first_or_diagnostic_visit_v131Jan24.rds")

crp_high_low_first_or_diagnostic_visit_v131Jan24$above_5 <- as.factor(crp_high_low_first_or_diagnostic_visit_v131Jan24$above_5)

df3 <- merge(crp_high_low_first_or_diagnostic_visit_v131Jan24, v4_clean_clinical_data_first_visit_15Dec23, by='id')


#now get graphs (for some reason not saved in earlier iteration.....)
#now compare clinical differences of interest between groups
clin123 = df3

#first plot clinical differences
labs <- c('<=5', '>5')
box_fills <- c('darkblue', 'darkred')

library(rstatix)

a1 <- ggplot(clin123, aes(x=above_5, y=hb_pawp_m)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('PAWP') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('Wedge pressure (mmHg)') +
  xlab('CRP group') +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)


a2 <- ggplot(clin123, aes(x=above_5, y=hb_rap_m)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('RAP') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('RAP (mmHg)') +
  xlab('CRP group') +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a3 <- ggplot(clin123, aes(x=above_5, y=bs_bmi)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('BMI') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('BMI (kg/m^2)') +
  xlab('CRP group') +
  # stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a4 <- ggplot(clin123, aes(x=above_5, y=egfr_mdrd)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('eGFR') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('eGFR (ml/min/1.73m^2)') +
  xlab('CRP group') +
  #  stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

clin1234 <- clin123 %>% filter(ep_1_type_6mwt == 'corridor')

a5 <- ggplot(clin1234, aes(x=above_5, y=ep_1_distance_meters)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('6MWD') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('6MWD (metres)') +
  xlab('CRP group') +
  # stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a6 <- ggplot(clin123, aes(x=above_5, y=lf_fvc_pc)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('FVC') +
  scale_fill_manual(values = box_fills, guide='none') +
  ylab('FVC (%pred)') +
  xlab('CRP group') +
  #stat_compare_means(method = "t.test") +
  scale_x_discrete(labels= labs) +
  theme_bw()+ 
  scale_y_continuous(expand=expansion(mult = c(0, .1)))
#  facet_wrap(~FinalPrimaryPHDiagnosis)

a7 <- ggplot(clin123, aes(x=above_5, y=hb_pvr_calc)) + 
  geom_boxplot(aes(fill = above_5), outlier.colour="black", outlier.shape=8, outlier.size=4) +
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

#=============================================



load("C:/Users/Gebruiker/Downloads/data_checked (4).RData")
sel <- data.checked %>% select(id, cbt_card_bnp_ngpl)
df3 <- merge(df3, sel, by='id')

wilcox.test(df3$bs_bmi ~ df3$above_5)
wilcox.test(df3$hb_pawp_m ~ df3$above_5)
wilcox.test(df3$hb_rap_m ~ df3$above_5)
wilcox.test(df3$lf_fvc_pc ~ df3$above_5)
wilcox.test(df3$cfe_rest_spo2 ~ df3$above_5)
wilcox.test(df3$egfr_mdrd ~ df3$above_5)


library(Publish)
library(broom)

CCI_score_per_patient <- readRDS("~/PhD/Projects/CRP ~ survival and BMI/CCI_score_per_patient.rds")

df3 <- merge(df3, CCI_score_per_patient, by='id')

df3$pvri <- (df3$hb_pap_m - df3$hb_pawp_m)/df3$hb_cardiac_index_value_1
#======================================================================================
#now get pvals required
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
output <- lapply(myvar_ch, function(x) do_chisq(df3, x, group_var))
combined_output1 <- bind_rows(output)


do_mw <- function(myvar) {
  mw_output <- wilcox.test(df3[[myvar]] ~ df3$above_5)
  mw_output <- (tidy(mw_output))
  mw_output$var <- myvar
  return(mw_output)
}

do_mw('age_diagnosis')

#now do for all
myvar_mw <- c('age_diagnosis', 'cbt_thyr_tsh_mupl', 'fev1_fvc', 'ep_1_distance_meters', 'cbt_thyr_freet4_pmolpl', 'bs_bmi', 'egfr_mdrd', 'cbt_haem_platelets_x10e9pl', 'cbt_haem_hb_gpl', 'cfe_rest_spo2', 'cfe_heart_rate', 'hb_pawp_m', 'hb_pap_d', 'CCI_adj_score', 'hb_pap_m', 'hb_pvr_calc', 'hb_rap_m', 'hb_cardiac_output_value_1', 'hb_cardiac_index_value_1', 'lf_fev1_pc', 'lf_fvc_pc', 'lf_kco_pc', 'cbt_card_ntprobnp_ngpl', 'cbt_card_bnp_ngpl', 'pvri')
myvar = myvar_mw
output <- lapply(myvar, do_mw)
combined_output2 <- bind_rows(output)

pvals <- rbind(combined_output1[,c(1,2,4,5)], combined_output2[,c(1:3,5)])

chisq.test(df3$diagnosis_verified, df3$above_5)
pvals$p.value <- ifelse(pvals$var == 'diagnosis_verified', 0.2634, pvals$p.value)

chisq.test(df3$hv_vasodilator_responder, df3$above_5)
pvals$p.value <- ifelse(pvals$var == 'hv_vasodilator_responder', 0.04977, pvals$p.value)

df4 <- df3 %>% filter(ep_1_type_6mwt == 'corridor')
wilcox.test(df4$ep_1_distance_meters ~ df4$above_5)

pvals$p.value <- data.table::fifelse(pvals$var == 'ep_1_distance_meters ', 1.289e-09, pvals$p.value)

#problem: pmerge 6MWD value still is wrong!
pvals$p.value2 <- ifelse(pvals$var == 'ep_1_distance_meters', 1.289e-09, NA)
pvals$p.value2 <- ifelse(is.na(pvals$p.value2), pvals$p.value, pvals$p.value2)


pvals$fdr <- p.adjust(pvals$p.value2, method='fdr')
pvals$fdradj <- round(pvals$fdr, digits = 10)

pmerge <- pvals %>% select(var, p.value2, fdr, fdradj)



a1 <- univariateTable(above_5 ~ diagnosis_verified + sex + age_diagnosis + bs_bmi + CCI_adj_score + functional_class  + cfe_rest_spo2 + cfe_heart_rate + ep_1_distance_meters + cbt_thyr_tsh_mupl + cbt_thyr_freet4_pmolpl + egfr_mdrd + cbt_haem_platelets_x10e9pl + cbt_haem_hb_gpl + cbt_card_ntprobnp_ngpl + cbt_card_bnp_ngpl + hb_pawp_m + hb_pap_d + hb_pap_m + hb_pvr_calc + pvri + hb_rap_m + hb_cardiac_output_value_1 + hb_cardiac_index_value_1 + hv_vasodilator_responder + lf_fev1_pc + lf_fvc_pc + fev1_fvc + lf_kco_pc, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = df3)
a2 <- summary(a1)

a2$order <- seq_len(nrow(a2))

a3 <- merge(a2, pmerge, by.x = 'Variable', by.y='var', all=T)
a3 <- a3[order(a3$order), ]


write.csv2(a3, 'differences_crp_groups_noClust.csv')


c1 <- univariateTable(above_5 ~ ep_1_distance_meters, summary.format = "median(x) [iqr(x)]",  column.percent = TRUE, compare.groups = TRUE, show.totals = TRUE, data = df4)
c2 <- summary(c1)



#compare comorbidity between groups
cleaned_comorbidity_dataV2_20Dec2023 <- readRDS("C:/Users/location/cleaned_comorbidity_dataV2_20Dec2023.rds")
m_hl <- crp_high_low_first_or_diagnostic_visit_v131Jan24

com <- merge(m_hl[,c(1,7)], cleaned_comorbidity_dataV2_20Dec2023, all=T)
com$above_5 <- as.factor(com$above_5)
com <- com %>% filter(!is.na(above_5))
com <- com[,1:3]

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

#now merge all of this by hand....
htn <- com %>% filter(htn == 'yes') %>% select(id, htn)
dm2 <- com %>% filter(dm2 == 'yes') %>% select(id, dm2)
hypothyr <- com %>% filter(hypothyr == 'yes') %>% select(id, hypothyr)
dyslip <- com %>% filter(dyslip == 'yes') %>% select(id, dyslip)
iheart <- com %>% filter(iheart == 'yes') %>% select(id, iheart)
copd <- com %>% filter(copd == 'yes') %>% select(id, copd)
asthma <- com %>% filter(asthma == 'yes') %>% select(id, asthma)
pe <- com %>% filter(pe == 'yes') %>% select(id, pe)
ckd <- com %>% filter(ckd == 'yes') %>% select(id, ckd)
obese <- com %>% filter(obese == 'yes') %>% select(id, obese)
mi <- com %>% filter(mi == 'yes') %>% select(id, mi)
osa <- com %>% filter(osa == 'yes') %>% select(id, osa)
iron <- com %>% filter(iron == 'yes') %>% select(id, iron)

df_list <- list(htn, dm2, hypothyr, dyslip, iheart, copd, asthma, pe, ckd, obese, mi, osa, iron)
df_try <- df_list %>% reduce(full_join, by='id')

com2 <- merge(com[,1:2], df_try, by='id', all=T)
com2[is.na(com2)] <- 'no'

com2[,3:15] <- lapply(com2[,3:15], as.factor)
com2 <- unique(com2)

library(Publish)
u1 <- univariateTable(above_5 ~ htn+ dm2+ hypothyr+ dyslip+ iheart+ copd+ asthma+ pe+ ckd+ obese+ mi+ osa+ iron, data=com2, column.percent = TRUE, show.totals = TRUE, compare.groups = TRUE)
u2 <- summary(u1)


#percentage_comrob
percom <- u2 %>% select(Variable, Level, `yes (n=381)`, `no (n=664)`)
percom <- percom %>% filter(Level == 'no')

percom[c('number_yes', 'percentage_yes')] <- str_split_fixed(percom$`yes (n=381)`, ' ', 2)
percom[c('number_no', 'percentage_no')] <- str_split_fixed(percom$`no (n=664)`, ' ', 2)

pernum <- percom %>% select(Variable, number_no, number_yes)
pernum$number_no <- as.numeric(pernum$number_no)
pernum$number_no_pos <- 664 - pernum$number_no
pernum$number_yes <- as.numeric(pernum$number_yes)
pernum$number_yes_pos <- 381 - pernum$number_yes

#also include a numeric BMI col
df3$obese <- ifelse(df3$bs_bmi >=30, 'obese by BMI', 'not obese')
df3$obese <- as.factor(df3$obese)

#get a column for OR
pernum$OR <- (pernum$number_yes_pos*pernum$number_no)/(pernum$number_yes*pernum$number_no_pos)
#calculate upper and lower CI: https://www.ncbi.nlm.nih.gov/books/NBK431098/
pernum$OR_CI_95_upper <- exp(log(pernum$OR) + 1.96*sqrt(((1/pernum$number_yes_pos) + (1/pernum$number_no) + (1/pernum$number_no_pos)) + (1/pernum$number_no)))
pernum$OR_CI_95_lower <- exp(log(pernum$OR) - 1.96*sqrt(((1/pernum$number_yes_pos) + (1/pernum$number_no) + (1/pernum$number_no_pos)) + (1/pernum$number_no)))


percom <- percom %>% select(Variable, percentage_no, percentage_yes)


percom <- percom %>% mutate(percentage_yes = str_replace_all(percentage_yes, "\\*|\\(|\\)", ""))
percom <- percom %>% mutate(percentage_no = str_replace_all(percentage_no, "\\*|\\(|\\)", ""))
percom$percentage_yes <- (100-as.numeric(percom$percentage_yes))
percom$percentage_no <- (100-as.numeric(percom$percentage_no))
rownames(percom) <- c('Hypertension', 'Type 2 Diabetes', 'Hypothyroidism', 'Dyslipidemia', 'Ischaemic heart disease', 'COPD', 'Asthma', 'Pulmonary embolism', 'CKD', 'Obesity', 'Myocardial infarction', 'OSA', 'Iron deficiency')



colnames(percom) <- c('var', 'Diagnostic CRP <=5mg/l', 'Diagnostic CRP >5mg/l')
#==============================================
#EdB continued here on 04-03-2024
#import smokingHx
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

smoke2 <- merge(smoke2, df[,c(1,315)], by.x='study_subject_oid', by.y='id')

#plot by group
library(Publish)
a1 <- univariateTable(above_5 ~ smoke_now, show.totals = TRUE, compare.groups = TRUE, data=smoke2)
a2 <- summary(a1)

a3 <- a2 %>% select(Level, `yes (n=372)`, `no (n=644)`)

a5 <- univariateTable(above_5 ~ obese, show.totals = TRUE, data = df3)
a6 <- summary(a5)

a3[c('number_yes', 'percentage_yes')] <- str_split_fixed(a3$`yes (n=372)`, ' ', 2)
a3[c('number_no', 'percentage_no')] <- str_split_fixed(a3$`no (n=644)`, ' ', 2)

a3 <- a3 %>% select(Level, percentage_no, percentage_yes)

a3 <- a3 %>% filter(Level %in% c('past_smoker', 'current_smoker'))

#give number the variables (to remove brackets etc.)
a3$percentage_no <- as.numeric(a3$percentage_no)
a3$percentage_yes <- as.numeric(a3$percentage_yes)
a3[1,2] <- 5.4
a3[2,2] <- 31.1
a3[1,3] <- 7.3
a3[2,3] <- 44.1

colnames(a3) <- c('var', 'Diagnostic CRP <=5mg/l', 'Diagnostic CRP >5mg/l')
rownames(a3) <- c('Current smoker', 'Past smoker')

permat2 <- rbind(percom, a3)
permat <- as.matrix(permat2[,2:3])

library(pheatmap)
pheatmap(permat, cluster_rows = FALSE, cluster_cols = FALSE, color=colorRampPalette(c('lightyellow', 'yellow', 'orange', 'darkorange', 'red', 'darkred'))(50))



pernum3 <- pernum %>% select(Variable, OR, OR_CI_95_upper, OR_CI_95_lower)
pernum3[14,c(1,2,3,4)] <- c('Obese by BMI', 2.7341, 3.5982, 2.0776)
pernum3[15,c(1,2,3,4)] <- c('Current smoker', 1.7150, 2.9464, 0.9982)
pernum3[16,c(1,2,3,4)] <- c('Past moker', 1.6848, 2.2153, 1.2821)
pernum3<- as.data.frame(pernum3)
write.csv2(pernum3, 'OR_comorbid_diseases_with_SE_crp.csv')

OR_comorbid_diseases_with_SE_crp <- read.csv("~/PhD/Projects/CRP ~ survival and BMI/OR_comorbid_diseases_with_SE_crp.csv", sep=";")

OR2 <- OR_comorbid_diseases_with_SE_crp

OR2$OR <- round(OR2$OR, digits = 4)
OR2$OR_CI_95_upper <- round(OR2$OR_CI_95_upper, digits = 4)
OR2$OR_CI_95_lower <- round(OR2$OR_CI_95_lower, digits = 4)

vars <- c('Hypertension', 'Type 2 diabetes', 'Hypothyrioidism', 'Dyslipidaemia', 'Ischaemic Heart Disease', 'COPD', 'Asthma', 'Pulmonary Embolism', 'CKD', 'Obese as comorbid diagnosis', 'Myocardial infarction', 'OSA', 'Iron deficiency', 'Obese by BMI', 'Current Smoker', 'Past smoker')

OR2$Variable <- vars

p <- 
  OR2 |>
  ggplot(aes(y = fct_rev(Variable))) + 
  theme_classic() +
  geom_point(aes(x=OR), shape=15, size=3) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Odds Ratio", y="") +
  geom_linerange(aes(xmin=OR_CI_95_lower, xmax=OR_CI_95_upper)) 

library(grid)
library(forestploter)

OR2$` ` <- paste(rep(" ", 20), collapse = "        ")


p <- forest(OR2[,c(2,3,6,4,5)],
            est = OR2$OR,
            lower = OR2$OR_CI_95_lower, 
            upper = OR2$OR_CI_95_upper,
            ci_column = 3,
            ref_line = 1,
            arrow_lab = c("Negative association", "Positive association"),
            xlim = c(-0.5, 4),
            ticks_at = c(0, 0.5, 1, 2, 3, 4),
            footnote = "Odds ratio for disease based on diagnostic CRP levels")

# Print plot
plot(p)
