#BMI change effect on 6MWD
#EdB started here on 16-02-2024
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
getwd()

#load packages
library(tidyverse)

#extract information from the latest datadump
load("interal_drive_location.RData")

#since CSCS is unreachable now
data_clean_cohort_used_for_crp <- readRDS("C:/Users/dropbox_data_location/data_clean_cohort_used_for_crp.rds")
data.clean.cohort <- data_clean_cohort_used_for_crp

#clean and prepare the DF
df <- data.clean.cohort %>% select(id, visit, ep_1_distance_meters, ep_1_type_6mwt, bs_bmi)
df2 <- df %>% select(id, visit, bs_bmi)

#code adapted from: https://stackoverflow.com/questions/68563299/calculate-percent-change-between-multiple-columns-of-a-data-frame
df2$visit <- as.integer(as.character((gsub('visit', '', df2$visit))))
df3 <- df2 %>% group_by(id) %>% arrange(visit) %>% mutate(percentual_change_BMI = (bs_bmi-lag(bs_bmi,1))/lag(bs_bmi,1)*100, change_bmi = bs_bmi-lag(bs_bmi,1))

#select labels for ones to plot
#first remove unlikely changes
df3 <- df3 %>% filter(percentual_change_BMI > -35 & percentual_change_BMI < 70)

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/Users/location/v4_clean_clinical_data_first_visit_15Dec23.rds")
labs <- v4_clean_clinical_data_first_visit_15Dec23 %>% select(id, age_diagnosis, diagnosis_verified)
#only adults but keep all causes of PAH in there
labs <- labs %>% filter(age_diagnosis >=18)

#merge with 6MWD % change
df4 <- df
df_walk <- df4 %>% filter(!is.na(ep_1_type_6mwt)) %>% group_by(id) %>% summarise(n_distinct(ep_1_type_6mwt))
#only take those patients for whom the walk test doesn't change
df_walk <- df_walk %>% filter(`n_distinct(ep_1_type_6mwt)` == 1)
keep_ids <- df_walk$id
df4 <- df4 %>% filter(id %in% keep_ids) %>% select(id, visit, ep_1_distance_meters) %>% group_by(id) %>% arrange(visit) %>% mutate(percentual_change_6MWD = (ep_1_distance_meters-lag(ep_1_distance_meters,1))/lag(ep_1_distance_meters,1)*100, abs_6mwd_change = ep_1_distance_meters-lag(ep_1_distance_meters,1))
#remove <0m walked (then test just not done it seems)
df4 <- df4 %>% filter(ep_1_distance_meters >=0)

#now merge
df3 <- merge(df3, df4, by=c('id', 'visit'))

#merge with labels of just adults
df3 <- merge(df3, labs, by='id')

#now plot
ggplot(df3, aes(x=percentual_change_BMI, y=percentual_change_6MWD)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change') +
  labs(y='Change in 6MWD (%)', x='Change in BMI (%)')


ggplot(df3, aes(x=change_bmi, y=abs_6mwd_change)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change') +
  labs(y='Change in 6MWD', x='Change in BMI')

#now compare 6MWD differences with stats
hist(df3$percentual_change_BMI)
hist(df3$percentual_change_6MWD)
#6MWD change is non-normal
hist(log(df3$percentual_change_6MWD))
#it is after log transforming
#now calculate pearson correlation
df3$log_mwd <- df3$percentual_change_6MWD + 100
df3$log_mwd <- log(df3$log_mwd)

cor.test(df3$percentual_change_BMI, df3$log_mwd, method="pearson", use="pairwise.complete.obs")
#this is NS

#Also test spearman
cor.test(df3$percentual_change_BMI, df3$percentual_change_6MWD, method="spearman", use="pairwise.complete.obs")
#this is NS

cor.test(df3$abs_6mwd_change, df3$change_bmi, method="spearman", use="pairwise.complete.obs")
#this is NS

#let's only take BMI changes of <20% to see what happens when we exclude severe outliers
df5 <- df3 %>% filter(percentual_change_BMI <20 & percentual_change_BMI >-20)

ggplot(df5, aes(x=percentual_change_BMI, y=percentual_change_6MWD)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change, outliers removed') +
  labs(y='Change in 6MWD (%)', x='Change in BMI (%)')

cor.test(df5$percentual_change_BMI, df5$percentual_change_6MWD, method="spearman", use="pairwise.complete.obs")
#this is still NS

#test significant 6MWD differences in population with BMI reduction over time (asssuming >5% is quite substantial)
df6 <- df3 %>% filter(percentual_change_BMI <=-5)

ggplot(df6, aes(x=percentual_change_BMI, y=percentual_change_6MWD)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change, outliers removed') +
  labs(y='Change in 6MWD (%)', x='Change in BMI (%)')

cor.test(df6$percentual_change_BMI, df6$percentual_change_6MWD, method="spearman", use="pairwise.complete.obs")
#this is sign


#=============================================================================
#repeat for absolute numbers
#test significant 6MWD differences in population with BMI reduction over time (asssuming >5% is quite substantial)
df6a <- df3 %>% filter(change_bmi <=-0.5)

ggplot(df6a, aes(x=change_bmi, y=abs_6mwd_change)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change') +
  labs(y='Change in 6MWD', x='Change in BMI')

cor.test(df6a$change_bmi, df6a$abs_6mwd_change, method="spearman", use="pairwise.complete.obs")
#this is almost sign
#=====================================================================

df7 <- df3 %>% filter(percentual_change_BMI >5)

ggplot(df7, aes(x=percentual_change_BMI, y=percentual_change_6MWD)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change, outliers removed') +
  labs(y='Change in 6MWD (%)', x='Change in BMI (%)')

cor.test(df7$percentual_change_BMI, df7$percentual_change_6MWD, method="spearman", use="pairwise.complete.obs")
#this is Nsign

#==================================================
#take positive BMI change
df7a <- df3 %>% filter(change_bmi >0.5)

ggplot(df7a, aes(x=change_bmi, y=abs_6mwd_change)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by BMI change, outliers removed') +
  labs(y='Change in 6MWD', x='Change in BMI')

cor.test(df7a$change_bmi, df7a$abs_6mwd_change, method="spearman", use="pairwise.complete.obs")
#this is sign
####################################################

test_df <- merge(data.clean.cohort, df6, by='id')
test_df$visit.x <- as.factor(test_df$visit.x)

kruskal.test(test_df$ep_1_distance_meters.x ~ test_df$visit.x)

ggplot(test_df, aes(x=visit.y, y=ep_1_distance_meters.x)) +
  geom_point() +
  stat_smooth() +
  ggtitle('6MWD change by visit in group with decreasing BMI') +
  labs(y='6MWD', x='Visit')
