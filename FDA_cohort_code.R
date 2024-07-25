# Libraries

library(dplyr)
library(survival)
library(survminer)
library(tidyverse)  # Pipe operator (%>%) and other commands
library(gtsummary)  
library(Publish)
library(broom)

# Load datasets
baseline_data <- read.csv("C:/Users/user.1/Desktop/baseline.csv")
followup_data <- read.csv("C:/Users/user.1/Desktop/followup.csv")
com <- read.csv("C:/Users/user.1/Desktop/comorbidities.csv")

#merging the FU data with comorbidities
updated <- merge(followup_data, com[-c(1)])

#adding needed variables to the FU dataset
FU <- updated %>%
  mutate(BMI = WEIGHT_KG/(HEIGHT_CM*0.01)^2) %>%
  mutate(BMI = round(BMI,1),
         survival_1y = if_else((DEATHDY >= 365), 365, as.double(DEATHDY)),
         survival_5y = if_else((DEATHDY >= 1825), 1825, as.double(DEATHDY)),
         YearsDeath = round(DEATHDY/365)) 

#adding BMI categories (from WHO)
FU <- FU %>%
  mutate(BMI_cat = case_when(BMI < 18.5 ~ "Underweight", #Underweight #1
                             BMI >= 18.5 & BMI <= 24.9  ~ "Normal weight", #Normal weight #2
                             BMI >= 25 & BMI <= 29.9 ~ "Pre-obesity", #Pre-obesity #3
                             BMI >= 30 & BMI <= 34.9 ~ "Obesity class-I", #Obesity class I #4
                             BMI >= 35 & BMI <= 39.9 ~ "Obesity class-II", #Obesity class II #5
                             BMI >= 40.0 ~ "Obesity class III")) #Obesity class III #6)

#checking the number of NAs
sum(is.na(FU$BMI)) #255
sum(is.na(FU$BMI_cat)) #255


# Recoding death
FU <- FU %>%
  mutate(status = if_else(DEATHFL == 'N', 0, #still alive
                          1)) #died
#reordening the BMI cat
new_order <- c("Obesity class III", "Obesity class-II", "Obesity class-I",
               "Pre-obesity", "Normal weight", "Underweight")

# Reorder the factor levels
FU$BMI_cat  <- factor(FU$BMI_cat, levels = new_order)
#checking the new order
table(FU$BMI_cat)

#apply the scaling
FU$pvr_scale <- scale(FU$PVR)
FU$rap_scale <- scale(FU$MRAP)
FU$age_scale <- scale(FU$AGE)

# Baseline data
baseline_data$DIAG_CD <- as.factor(baseline_data$DIAG_CD)
baseline_data$FC <- as.factor(baseline_data$FC)

#adding needed variables
bsl <- baseline_data %>%
  mutate(BMI = WEIGHT_KG/(HEIGHT_CM*0.01)^2,
         survival_1y = if_else((DEATHDY >= 365), 365, as.double(DEATHDY)),
         survival_5y = if_else((DEATHDY >= 1825), 1825, as.double(DEATHDY)),
         YearsDeath = round(DEATHDY/365),
         SEX = factor(SEX),
         PVRI = ((MPAP - MPCWP)/ CI),
         TPR = MPAP/CO)

#round
bsl$BMI <- round(bsl$BMI,1)

# BMI categories for baseline data
bsl <- bsl %>%
  mutate(BMI_cat = factor(case_when(BMI < 18.5 ~ "Underweight", #Underweight #1
                                    BMI >= 18.5 & BMI <= 24.9  ~ "Normal weight", #Normal weight #2
                                    BMI >= 25.0 & BMI <= 29.9 ~ "Pre-obesity", #Pre-obesity #3
                                    BMI >= 30.0 & BMI <= 34.9 ~ "Obesity class-I", #Obesity class I #4
                                    BMI >= 35.0 & BMI <= 39.9 ~ "Obesity class-II", #Obesity class II #5
                                    BMI >= 40.0 ~ "Obesity class III"),#Obesity class III #6
  ))

#checking distribution of baseline data
#selecting specific numeric variables
numericsTable <- bsl[,c(4,76,21,22,24,19,20,10)]

data_long <- numericsTable %>%
  pivot_longer(colnames(numericsTable)) %>%
  as.data.frame()

ggp1 <- ggplot(data_long, aes(x = value)) +
  geom_histogram(aes(y = ..density..)) +
  geom_density(col = "#1b98e0", size = 1)+
  facet_wrap(~name, scales = "free")
ggp1

#checking normality with Shapiro
#Shapiro-Wilk test
#H0 : the variable follow a normal distribution
#H1 : the variable does NOT follow a normal distribution
apply(numericsTable,2,shapiro.test)

#checking QQ plots
apply(numericsTable,2,ggqqplot)

### the p-values < 0.05 are implying that the distribution of the data are significantly different from normal distribution. In other words, we can't assume the normality.Visuals agreed

#excluding NAs from BMI categories
BMI_no_NA <- bsl %>%
  drop_na(BMI_cat)


# Table with baseline characteristics
uv1 <- univariateTable(BMI_cat ~ SEX + AGE + DIAG_CD + PVR + PVRI + TPR + CI + 
                         CO + MRAP + MPAP + MPCWP +  EGFR + SIXMWT_D + FC, data=BMI_no_NA, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE, summary.format = 'median(x) [iqr(x)]')#we dont have FVCp +
uv2 <- summary(uv1) #, show.missing = c("never")
uv2

# Saving to pc
write.csv(uv2,"C:/Users/user1/bsl_table.csv")

#table baseline differences by diagnosis
a1 <-univariateTable(DIAG_CD ~ SEX + AGE + BMI +  MRAP + MPAP + MPCWP + CO + PVR + SIXMWT_D + FC, data=bsl, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
a2 <- summary(a1, show.missing = c("never"))

# Eck krustal function
do_kruskal2 = function(myvar) {
  output_krusk <- kruskal.test(BMI_no_NA[[myvar]] ~ BMI_no_NA$BMI_cat)
  output_krusk = tidy(output_krusk)
  output_krusk$var <- myvar
  return(output_krusk)
}


do_kruskal2('AGE')

#this works
myvar = c('SEX','AGE','DIAG_CD','PVR','PVRI','TPR','CI','CO','MRAP','MPAP','MPCWP','EGFR','SIXMWT_D')
o1 <- lapply(myvar, do_kruskal2)
o2 <- bind_rows(o1)
o2

# assess treatment effects
sixbsl <- bsl %>%
  select(USUBJID,BMI_cat, BMI, SIXMWT_D, DIAG_CD, AGE)

sixfu <- FU %>%
  select(USUBJID, SIXMWT_D,BMI_cat, AGE) %>%
  rename(IDFU = USUBJID,
         BMI_catFU = BMI_cat,
         SIXFU = SIXMWT_D,
         AGEFU = AGE)

finsix <- cbind(sixbsl,sixfu)

#calculate the 6mwt difference between FU and Bsl
finsix$Sixdifference <- finsix$SIXFU - finsix$SIXMWT_D 

# Plotting 6mwt difference
# Histogram with density plot
ggplot(finsix, aes(x=Sixdifference)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

#reordering the BMI cat levels
new_order <- c("Obesity class III", "Obesity class-II", "Obesity class-I",
               "Pre-obesity", "Normal weight", "Underweight")

# Reorder the factor levels
finsix$BMI_cat  <- factor(finsix$BMI_cat, levels = new_order)
table(finsix$BMI_cat)

#plotting difference in 6mwd between BMI categories
h <- finsix %>%
  drop_na(BMI_cat) %>%
  ggplot(aes(x=BMI_cat, y=Sixdifference)) + 
  geom_boxplot(aes(fill = BMI_cat), outlier.colour="black", outlier.shape=8, outlier.size=4) +
  ggtitle('Difference in 6MWD improvement between BMI groups') +
  ylab('6MWD improvement (m)') +
  xlab('BMI group') +
  theme_bw() +
  stat_compare_means(method = "anova")

h <- h + guides(fill=guide_legend(title="Weight"))
h

# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("SIXMWD_diff_6.15.24.png", h, height = 6, width = 8)

#COX 
FU$BMI_cat <- relevel(FU$BMI_cat, ref = "Normal weight")

#excluding missing data
FU_BMI_no_NA <- FU %>%
  drop_na(BMI_cat)

soa <- Surv(FU_BMI_no_NA$DEATHDY/365, FU_BMI_no_NA$status == 1, type = 'right')
sob <- survfit(soa ~BMI_cat, data = FU_BMI_no_NA)

# CoxPH
cox_bmi6 <- coxph(soa ~ BMI_cat + SEX + age_scale + pvr_scale  + rap_scale , data=FU_BMI_no_NA)
ggforest(cox_bmi6, data=FU_BMI_no_NA)
summary_cox3 <- summary(cox_bmi6)

#Output table
model1 <- tbl_regression(cox_bmi6, intercept = TRUE,  exp = TRUE)
as_gt(model1) %>% 
  gt::tab_header("Cox proportional hazards regression model") %>% 
  gt::tab_options(table.align='left')


#KM curve
p <- ggsurvplot(sob, data=FU_BMI_no_NA, 
                pval=TRUE, 
                pval.method = TRUE, 
                risk.table = TRUE, 
                conf.int = TRUE, 
                xlab='Time- (years)', 
                title='Survival differences based on BMI', 
                xlim=c(0,2), 
                ylim=c(0.6,1), 
                legend.title='Group', 
                break.x.by=1, 
                legend.labs=c("obesity class-III", "obesity class-II", "obesity class-I", "pre-obesity", "normal weight", "underweight"))
p
