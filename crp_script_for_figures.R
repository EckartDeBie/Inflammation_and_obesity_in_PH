#figures for CRP
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


log_transformed_crp_levels_clusters <- readRDS("C:/Users/Gebruiker/Documents/PhD/Projects/CRP ~ survival and BMI/log_transformed_crp_levels_clusters.rds")

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")

clin1 <- merge(log_transformed_crp_levels_clusters, v4_clean_clinical_data_first_visit_15Dec23, by.x='name', by.y='id')
clin1$value <- as.factor(clin1$value)

labs <- c('Normal CRP', 'High CRP')

p1 <- ggplot(clin1, aes(x=value, y=hb_pawp_m, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Wedge pressure') +
  ylab('Wedge pressure (mmHg)') +
  xlab('cluster') +
  theme_bw() +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs)

p2 <- ggplot(clin1, aes(x=value, y=hb_rap_m, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('RAP') +
  ylab('RAP (mmHg)') +
  xlab('cluster') +
  theme_bw() +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs)

p3 <- ggplot(clin1, aes(x=value, y=bs_bmi, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('BMI') +
  ylab('BMI (kg/m^2)') +
  xlab('cluster') +
  theme_bw() +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs)

p4 <- ggplot(clin1, aes(x=value, y=egfr_mdrd, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('eGFR') +
  ylab('eGFR (MDRD;ml/min/1.73m^2)') +
  xlab('cluster') +
  theme_bw() +
  #stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels= labs)

b <- clin1 %>% filter(ep_1_type_6mwt == 'corridor')

p5 <- ggplot(b, aes(x=value, y=ep_1_distance_meters, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('6MWD') +
  ylab('6MWD (m)') +
  xlab('cluster') +
  #stat_compare_means(method = "wilcox.test") +
  theme_bw() +
  scale_x_discrete(labels= labs)


p6 <- ggplot(clin1, aes(x=value, y=lf_fvc_pc, fill=value)) + 
  geom_boxplot(fill = c('darkblue', 'darkred'), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('FVC') +
  ylab('FVC (%pred)') +
  xlab('cluster') +
 # stat_compare_means(method = "wilcox.test") +
  theme_bw() +
  scale_x_discrete(labels= labs)


library(patchwork)
p1 + p2 + p3 +p4 + p5 +p6 + plot_layout(ncol=6)



#compare AAB levels!
#import the AABs
log.transformed.autoAbs.of.ctrls.and.pts.v1.3March2021 <- read.csv2("~/Not PhD/My publications/Data new analyses March 2021/log-transformed autoAbs of ctrls and pts v1 3March2021.csv")
autoimmunity_dataframe <- read_excel("~/Not PhD/TBR PAH cohort study/Origional data with minor transformations/Original data/autoimmunity.dataframe.xlsx")
autoimmunity_dataframe <- as.data.frame(autoimmunity_dataframe)
colnames(autoimmunity_dataframe) <- autoimmunity_dataframe[1,]
autoimmunity_dataframe <- autoimmunity_dataframe[-1,]
ai <- autoimmunity_dataframe %>% select(id_oc, id_cohort.x)
abs <- merge(log.transformed.autoAbs.of.ctrls.and.pts.v1.3March2021, ai, by.x='sample', by.y='id_cohort.x')

ab2 <- merge(abs, clin1[,c(1,2)], by.x = 'id_oc', by.y='name')
ab2 <- unique(ab2)

#retains 171 pts

names(ab2) <- c('id', 'sample', 'x', 'Cardiolipin', 'CENP-B', 'H2a/F2a & H4a/F2a1', 'Histone IIA', 'Jo-1', 'La/SS-B', 'Mi-2b', 'MPO', 'Proteinase-3', 'PDH', 'RNP complex', 'Ro/SS-A', 'Scl-34', 'Scl-70', 'Smith', 'Thyroglobulin', 'TPO', 'Transglutaminase', 'U1-snRNP68', 'sex', 's', 'age', 'pah', 'above_5')
#now group AABs
ab2 <- ab2 %>% pivot_longer(4:22, names_to='Autoantibody', values_to = 'AAB_level')

ab2$above_5 <- ifelse(ab2$above_5 == '1', 'High CRP', 'Normal CRP')

library(ggpubr)
levels(ab2$above_5)
ab2$above_5 <- as.factor(ab2$above_5)
ab2$above_5 <- relevel(ab2$above_5, ref=2)

ggplot(ab2, aes(x = above_5, y = AAB_level, fill = above_5)) + 
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 4) +
  ggtitle('Difference in Autoantibody level between CRP clusters') +
  ylab('Log-normalised AAB level') +
  xlab('CRP cluster') +
  theme_bw() +
  #stat_compare_means(method = 'wilcoxon.test', label = "p.format", vjust = 1) +
  facet_wrap(~ Autoantibody, scales = "free_y") +
  scale_fill_manual(values = c('darkblue', 'darkred'))

results <- ab2 %>%
  group_by(Autoantibody) %>%
  do(tidy(wilcox.test(AAB_level ~ above_5, data = .)))

# Print the results
print(results)

