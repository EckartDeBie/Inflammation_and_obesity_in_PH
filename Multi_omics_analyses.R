#MultiOmics for paper
#RNA-seq data
getwd()
library(tidyverse)
library(readxl)
#now load in the relevant files
load('/filepath/res.tpm.RData')
#transpose
tmp_trans <- t(res.tpm)
tmp_trans <- as.data.frame(tmp_trans)

#load the groups
crp_groups <- readRDS('/filepath/crp_high_low_first_or_diagnostic_visit_v131Jan24.rds')
library(readxl)
id_conv <- read_excel('/filepath/RNA_ID_convert.xlsx')
id_key <- id_conv

View(id_conv)
id_conv <- id_conv %>% select(SampleID, OpenClinicaID)
id2 <- id_conv
id_conv <- merge(id_conv, crp_groups, by.x='OpenClinicaID', by.y='id')
View(id_conv)
#312 patients overlap with RNA-seq data and CRP high vs low groups

#merge this with the RNA-Seq data
tmp_trans$SampleID <- rownames(tmp_trans)
#now just select the genes of interest

#first select genes of interest:
genes <- read_excel('/filepath/Genes_of_interest_inflammation.xlsx')
interst <- c(genes$gene_name)

#do analysis
names.use <- names(tmp_trans)[names(tmp_trans) %in% interst]
tmp_try <- tmp_trans[, names.use]
tmp_try$SampleID <- rownames(tmp_try)

tmp_try <- merge(tmp_try, id_conv, by='SampleID')
View(tmp_try)

#now see if there are significant differences
#first plot histograms
#310 patients left with RNA-Seq and labels
#extract the file for now 
saveRDS(tmp_try, file='RNA_expression_of_inflammatory_genes.rds')


#=====================================================================================
#re-run:
#re-do joint model
#remove duplicate OC IDs
tmp_try <- tmp_try %>% filter(OpenClinicaID != "") %>% distinct(OpenClinicaID, .keep_all = TRUE)    
#310 retained


clin_dat_at_diagnosis_clean <- readRDS("~/filepath/clin_dat_at_diagnosis_clean.rds")
follow_up_data <- readRDS("~/filepath/follow_up_data.rds")
#TOTAL MODEL OF METABOLITES
clina <- clin_dat_at_diagnosis_clean %>% select(id, bs_bmi, cbt_inflammation_crp_mgpl)
clina <- na.omit(clina)

clinb <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex)
clinb <- na.omit(clinb)

follow_up_data <- readRDS("~/filepath/follow_up_data.rds")
fu <- follow_up_data %>% select(id, id_cohort) %>% unique()
clina <- merge(clina, fu, by='id')
clina <- na.omit(clina)
clina$log_crp <- log(clina$cbt_inflammation_crp_mgpl)
clina <- merge(clina, clinb, by='id')


tmp_try <- tmp_try[,-c(88:93,1)]
p2 <- merge(clina, tmp_try, by.x='id', by.y='OpenClinicaID')
#256 patients retained

numeric_cols <- p2 %>% select(-id, -bs_bmi, -cbt_inflammation_crp_mgpl, -log_crp, - id_cohort, -age_diagnosis, -sex) %>% select(where(is.numeric))

library(broom)
results <- data.frame()

#get loop
for (col_name in colnames(numeric_cols)) {
  model <- lm(log(p2[[col_name]] + 10) ~ bs_bmi + log_crp + age_diagnosis + sex, data = p2)
  model_summary <- summary(model)$coefficients
  log_crp_coef <- model_summary["log_crp", "Estimate"]
  log_crp_pval <- model_summary["log_crp", "Pr(>|t|)"]
  bs_bmi_coef <- model_summary["bs_bmi", "Estimate"]
  bs_bmi_pval <- model_summary["bs_bmi", "Pr(>|t|)"]
  results <- rbind(results, data.frame(
    variable = col_name,
    log_crp_coef = log_crp_coef,
    log_crp_pval = log_crp_pval,
    bs_bmi_coef = bs_bmi_coef,
    bs_bmi_pval = bs_bmi_pval
  ))
}

results$fdr_crp <- p.adjust(results$log_crp_pval, method='fdr')
results$fdr_bmi <- p.adjust(results$bs_bmi_pval, method='fdr')

results$log10crp <- -log10(results$log_crp_pval)
results$log10bm <- -log10(results$bs_bmi_pval)

getwd()
write.csv2(results, file='transcriptomics_linear_model_crp_bmi.csv')

#now plot!
ggplot(results, aes(x=log10bm, y=log10crp)) +
  geom_point(aes(color = log10bm > 2.6 | log10crp > 4), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 10.5) +
  ylim(0, 10.5) +
  ggtitle('Differences in proteomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 4, linetype='dotted', col='red') +
  geom_vline(xintercept = 2.6, linetype='dotted', col='red') +
  annotate("text", x = 3.3, y = 4.1, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(results, (log10bm > 2.6 | log10crp > 4)), 
                  aes(x=log10bm, y=log10crp, label=variable), size = 3.5) +
  labs(color = "FDR significant")

#======================================
#use groups
weight <- readRDS('/filepath/BMI_category_per_patient.rds')
clina <- merge(clina, weight, by='id')
clina$above_5 <- ifelse(clina$cbt_inflammation_crp_mgpl >5, 'yes', 'no')

p3 <- merge(clina, tmp_try[,-c(96:101)], by.x='id', by.y='OpenClinicaID')
p3$above_5 <- as.factor(p3$above_5)

numeric_cols <- p3 %>% select(-id, -bs_bmi, -cbt_inflammation_crp_mgpl, -log_crp, - id_cohort, -age_diagnosis, -sex, -weight, -above_5) %>% select(where(is.numeric))
results <- data.frame()

# do loop


# Loop over each column in numeric_cols
for (col_name in colnames(numeric_cols)) {
  # Fit the linear model with log transformation for the response variable
  model <- lm(log(p3[[col_name]] + 10) ~  above_5, data = p3)
  
  # Extract coefficients and p-values from the model summary
  model_summary <- summary(model)$coefficients
  
  # Initialize coefficients and p-values for weight and above_5
  #weight_coefs <- model_summary[grep("^weight", rownames(model_summary)), "Estimate"]
  #weight_pvals <- model_summary[grep("^weight", rownames(model_summary)), "Pr(>|t|)"]
  
  # Make sure weight_coef and weight_pval have the same number of rows, and fill with NAs if needed
  # weight_coefs <- c(weight_coefs, rep(NA, max(0, length(levels(p3$weight)) - length(weight_coefs))))
  # weight_pvals <- c(weight_pvals, rep(NA, max(0, length(levels(p3$weight)) - length(weight_pvals))))
  
  # Loop through each level of 'above_5' (excluding the reference level) and extract coefficients and p-values
  for (level in levels(p3$above_5)[-1]) {  # Skip the reference level
    above_5_level <- paste0("above_5", level)  # Construct the variable name for each level
    
    # Check if this level exists in model_summary
    above_5_coef <- if (above_5_level %in% rownames(model_summary)) {
      model_summary[above_5_level, "Estimate"]
    } else {
      NA
    }
    above_5_pval <- if (above_5_level %in% rownames(model_summary)) {
      model_summary[above_5_level, "Pr(>|t|)"]
    } else {
      NA
    }
    
    # Create a row of results for the current variable and `above_5` level
    results_row <- data.frame(
      variable = col_name,
      above_5_level = level,
      above_5_coef = above_5_coef,
      above_5_pval = above_5_pval,
      #weight_coef = paste(weight_coefs, collapse = ", "),  # Combine weight coefficients as a comma-separated string
      #weight_pval = paste(weight_pvals, collapse = ", ")   # Combine weight p-values as a comma-separated string
    )
    
    # Append the results row to the dataframe
    results <- rbind(results, results_row)
  }
}


results$fdr_crp <- p.adjust(results$above_5_pval, method='fdr')
results$fdr_bmi <- p.adjust(results$bs_bmi_pval, method='fdr')

#=========================================================================================
# Initialize an empty dataframe to store results
# Initialize an empty dataframe with the correct column names
results <- data.frame(
  variable = character(),
  above_5_level = character(),
  above_5_coef = numeric(),
  above_5_pval = numeric(),
  intercept_coef = numeric(),
  intercept_pval = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each column in numeric_cols
for (col_name in colnames(numeric_cols)) {
  # Fit the linear model with log transformation for the response variable
  model <- lm(log(p3[[col_name]] + 10) ~ above_5, data = p3)
  
  # Extract coefficients and p-values from the model summary
  model_summary <- summary(model)$coefficients
  
  # Extract the intercept (reference level) coefficient and p-value
  intercept_coef <- model_summary["(Intercept)", "Estimate"]
  intercept_pval <- model_summary["(Intercept)", "Pr(>|t|)"]
  
  # Loop through each level of 'above_5' (excluding the reference level) and extract coefficients and p-values
  for (level in levels(p3$above_5)[-1]) {  # Skip the reference level
    above_5_level <- paste0("above_5", level)  # Construct the variable name for each level
    
    # Check if this level exists in model_summary
    above_5_coef <- if (above_5_level %in% rownames(model_summary)) {
      model_summary[above_5_level, "Estimate"]
    } else {
      NA
    }
    above_5_pval <- if (above_5_level %in% rownames(model_summary)) {
      model_summary[above_5_level, "Pr(>|t|)"]
    } else {
      NA
    }
    
    # Create a row of results for the current variable and `above_5` level
    results_row <- data.frame(
      variable = col_name,
      above_5_level = level,
      above_5_coef = above_5_coef,
      above_5_pval = above_5_pval,
      intercept_coef = intercept_coef,  # Include the reference level (intercept) coefficient
      intercept_pval = intercept_pval,   # Include the reference level (intercept) p-value
      stringsAsFactors = FALSE
    )
    
    # Append the results row to the dataframe
    results <- rbind(results, results_row)
  }
}

results$crp_fdr <- p.adjust(results$above_5_pval, method='fdr')


#proteomics data extraction for time series clustering
library(tidyverse)
library(ggpubr)

setwd('/home/emddd2/projects/PHomics')
getwd()

#import metabolomics
metablomics_normalisd_z_score_metabolomics <- readRDS("~/projects/Data_for_cynapse/Omics/Metabolomics/metablomics_normalisd_z_score_metabolomics.rds")
weight <- readRDS('/filepath/BMI_category_per_patient.rds')
weight <- weight %>% filter(!is.na(weight))
follow_up_data <- readRDS("~/filepath/follow_up_data.rds")

weight <- merge(weight, follow_up_data[,c(1,6)], by='id')
weight <- weight %>% filter(!is.na(id_cohort))
weight <- unique(weight)
#720 patients left
#list from: https://pmc.ncbi.nlm.nih.gov/articles/PMC5311903/
mito_metab <- c('citrate', 'cis.aconitate', 'fumarate', 'malate', 'pyruvate', 'succinate', 'glutamate', 'N6.succinyladenosine', 'succinylcarnitine', 'acetylcarnitine', 'id', 'weight')

metab <- merge(weight, metablomics_normalisd_z_score_metabolomics, by.x='id_cohort', by.y='SUBJECT_ID')
#238 patients with both
metab <- metab %>% select(all_of(mito_metab))
metab <- metab %>% filter(!is.na(weight))
#230 patients left

m2 <- metab %>% pivot_longer(1:10, names_to = 'metabolite', values_to = 'level')
m2 <- m2 %>% filter(!is.na(weight))
hist(m2$level)
#is normal!

m2$weight <- factor(m2$weight, levels = c("underweight", "normal weight", "pre-obesity", 
                                          "obesity class-I", "obesity class-II", "obesity class-III"), 
                    ordered = TRUE)

#now plot
ggplot(m2, aes(x = weight, y = level, fill=weight)) +
  geom_boxplot() +
  facet_wrap(~ metabolite, scales = 'free') +
  stat_compare_means(method = "anova", label = "p.format") +
  theme_minimal() +
  labs(y = "metabolite level", x = "weight group", title = "metabolite levels")


#################################################################################################################
#get all differences
metabol <- merge(weight, metablomics_normalisd_z_score_metabolomics, by.x='id_cohort', by.y='SUBJECT_ID')
#order the factor
metabol$weight <- factor(metabol$weight, levels = c("underweight", "normal weight", "pre-obesity", 
                                                    "obesity class-I", "obesity class-II", "obesity class-III"), 
                         ordered = TRUE)

numeric_cols <- metabol[, sapply(metabol[, 4:1419], is.numeric)]
numeric_cols <- numeric_cols %>% select(id, weight, where(is.numeric))

#get distribution
num <- numeric_cols%>% pivot_longer(3:1311, names_to = 'met', values_to = 'val')
hist(num$val, breaks=250)
#roughly normal distribution, slight skew but large numbers


numeric_cols$weight <- factor(numeric_cols$weight, levels = c("underweight", "normal weight", "pre-obesity", 
                                                              "obesity class-I", "obesity class-II", "obesity class-III"), 
                              ordered = TRUE)

numeric_cols <- metabol %>% select(-id, -weight) %>% select(where(is.numeric))

#=============================================================================
#add in CRP and BMI
clin_dat_at_diagnosis_clean <- readRDS("~/filepath/clin_dat_at_diagnosis_clean.rds")
clin1 <- clin_dat_at_diagnosis_clean %>% select(id, cbt_inflammation_crp_mgpl, bs_bmi)

hist(clin1$cbt_inflammation_crp_mgpl)
#is not normal
clin1$log_crp <- log(clin1$cbt_inflammation_crp_mgpl)
hist(clin1$log_crp)
#much better and normal now
hist(clin1$bs_bmi)
#is normal


#now get residuals of metabolomics for each CRP
clin1 <- clin1 %>% filter(!is.na(log_crp))
metabol2 <- merge(metabol, clin1, by='id')
#202 retained if CRP present for adjustment

#now get residuals
# Initialize an empty dataframe to store residuals
residuals_df <- data.frame(id = metabol2$id)

# Calculate residuals for each numeric column in `numeric_cols`
for (col_name in colnames(numeric_cols)) {
  # Fit the linear model with `log_crp` and `bs_bmi` as predictors
  model <- lm(metabol2[[col_name]] ~ log_crp, data = metabol2)
  
  # Extract residuals and add them as a new column in `residuals_df`
  residuals_df[[col_name]] <- residuals(model)
}

#now merge these with weights groups
residuals_df <- merge(metabol2[,c(1,3)], residuals_df, by='id')

#metabol <- residuals_df
#======================================================================

#Perform ANOVA for each numeric column and extract coefficients
anova_results <- lapply(names(numeric_cols), function(col_name) {
  formula <- as.formula(paste(col_name, "~ weight"))
  aov_result <- aov(formula, data = metabol)
  
  # Extract p-value
  p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
  
  # Extract coefficients
  coefficients <- coef(aov_result)
  
  # Combine into a data frame
  data.frame(
    Variable = col_name,
    Intercept = coefficients[1],
    Slope = coefficients[2],
    p_value = p_value
  )
})

# Combine results into a single dataframe and include FDR
anova_summary_df <- do.call(rbind, anova_results)
anova_summary_df$FDR <- p.adjust(anova_summary_df$p_value, method='fdr')
a2 <- anova_summary_df %>% filter(FDR <0.05)
#===========================================================================
#now do CRP
crp_groups <- readRDS('/filepath/CRP analyses/crp_high_low_first_or_diagnostic_visit_v131Jan24.rds')
mito_metab2 <- c('citrate', 'cis.aconitate', 'fumarate', 'malate', 'pyruvate', 'succinate', 'glutamate', 'N6.succinyladenosine', 'succinylcarnitine', 'acetylcarnitine', 'id', 'above_5')


metab2 <- merge(crp_groups, metablomics_normalisd_z_score_metabolomics, by.x='id_cohort', by.y='SUBJECT_ID')
#204 patients with both
metab2 <- metab2 %>% select(all_of(mito_metab2))
metab2 <- metab2 %>% filter(!is.na(above_5))
#204 patients left

m2 <- metab2 %>% pivot_longer(1:10, names_to = 'metabolite', values_to = 'level')

#now plot
ggplot(m2, aes(x = above_5, y = level, fill=above_5)) +
  geom_boxplot() +
  facet_wrap(~ metabolite, scales = 'free') +
  stat_compare_means(method = "t.test", label = "p.format") +
  theme_minimal() +
  labs(y = "metabolite level", x = "CRP group", title = "metabolite levels")


#=======================================================================================
#also total analysis
#get all differences
metabol <- merge(crp_groups, metablomics_normalisd_z_score_metabolomics, by.x='id_cohort', by.y='SUBJECT_ID')
metabol <- metabol[,-c(1,3:6)]
#order the factor
numeric_cols <- metabol %>% select(-id, -above_5) %>% select(where(is.numeric))

# Ensure `above_5` is a factor
metabol$above_5 <- as.factor(metabol$above_5)
###################################################################################
#include adjustment
clin1 <- clin_dat_at_diagnosis_clean %>% select(id, cbt_inflammation_crp_mgpl, bs_bmi)

#now get residuals of metabolomics for each CRP
clin1 <- clin1 %>% filter(!is.na(bs_bmi))
metabol2 <- merge(metabol, clin1, by='id')
#197 retained if CRP present for adjustment

#now get residuals
# Initialize an empty dataframe to store residuals
residuals_df <- data.frame(id = metabol2$id)

# Calculate residuals for each numeric column in `numeric_cols`
for (col_name in colnames(numeric_cols)) {
  # Fit the linear model with `log_crp` and `bs_bmi` as predictors
  model <- lm(metabol2[[col_name]] ~ bs_bmi, data = metabol2)
  
  # Extract residuals and add them as a new column in `residuals_df`
  residuals_df[[col_name]] <- residuals(model)
}

#now merge these with weights groups
residuals_df <- merge(metabol2[,c(1,2)], residuals_df, by='id')

#metabol <- residuals_df
##############################################################################

#do test
t_test_results <- lapply(names(numeric_cols), function(col_name) {
  # Check if both groups in 'above_5' have more than one unique value in the numeric column
  if(length(unique(metabol[[col_name]][metabol$above_5 == unique(metabol$above_5)[1]])) > 1 && 
     length(unique(metabol[[col_name]][metabol$above_5 == unique(metabol$above_5)[2]])) > 1) {
    
    # Create the formula for the t-test
    formula <- as.formula(paste(col_name, "~ above_5"))
    
    # Perform t-test
    t_test_result <- t.test(formula, data = metabol)
    
    # Extract p-value
    p_value <- t_test_result$p.value
    
    # Calculate the slope (mean difference) between groups
    group_means <- tapply(metabol[[col_name]], metabol$above_5, mean)
    slope <- diff(group_means) # Difference in means between the two groups
    
    # Return results as a data frame
    return(data.frame(Variable = col_name, p_value = p_value, slope = slope))
  } else {
    # If there's no variability, return NA for p-value and slope
    return(data.frame(Variable = col_name, p_value = NA, slope = NA))
  }
})

# Combine results into a single data frame
t_test_results <- do.call(rbind, t_test_results)

# Combine results into a single dataframe
t_test_summary_df <- t_test_results
# Adjust p-values using FDR (False Discovery Rate) correction
t_test_summary_df$FDR <- p.adjust(t_test_summary_df$p_value, method='fdr')

# Filter for significant results (FDR < 0.05)
significant_t_tests <- t_test_summary_df %>% filter(FDR < 0.05)

#test overlap
abc123 <- merge(a2, significant_t_tests, by='Variable')
#32 overlap ~50% (0% overlap if adjusted)

#now plot
t1 <- t_test_summary_df 
t1$logpval_crp <- -log10(t1$p_value)

t2 <- anova_summary_df
t2$logval_bmi <- -log10(t2$p_value)

t3 <- merge(t1, t2, by='Variable')

plot(t3$logpval_crp, t3$logval_bmi)

#========================================================
#get plot for Chris
ta <- t_test_summary_df 
names(ta) <- c('Variable', 'CRP_pval', 'CRP_slope', 'CRP_FDR')

tb <- anova_summary_df
names(tb) <- c('Variable', 'bmi_intercept', 'BMI_slope', 'BMI_pval', 'BMI_FDR')

tab <- merge(ta, tb, by='Variable')
tab$significance <- ifelse(tab$BMI_FDR < 0.05 & tab$CRP_FDR <0.05, 'both', NA)
tab$significance <- ifelse(tab$BMI_FDR < 0.05 & tab$CRP_FDR >=0.05, 'BMI only', tab$significance)
tab$significance <- ifelse(tab$BMI_FDR >=0.05 & tab$CRP_FDR <0.05, 'CRP only', tab$significance)
tab$significance <- ifelse(tab$BMI_FDR >=0.05 & tab$CRP_FDR >=0.05, 'none', tab$significance)

tab <- tab %>% filter(!is.na(significance))

ggplot(tab, aes(x=BMI_slope, y=CRP_slope, colour=significance)) +
  geom_point(size = 2, shape = 19) +  
  theme_bw() +
  xlim(-1.5, 2) +
  ylim(-1.5, 2) +
  ggtitle('Differences in metabolomic signature') +
  xlab('Slope for BMI groups') +
  ylab('Slope for CRP groups') +
  labs(color = "FDR significance") +
  scale_color_manual(values = c("CRP only" = "red", "none" = "black", "BMI only" = "blue", "both" = "purple"))

#================================================
#now plot
library(ggrepel)

t3 <- t3 %>% filter(!is.na(logpval_crp))

plot_scaled <- ggplot(t3, aes(x=logval_bmi, y=logpval_crp)) +
  geom_point(aes(color = logval_bmi > 3.3 | logpval_crp > 3.3), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 8) +
  ylim(0, 8) +
  ggtitle('Differences in metabolomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 3.3, linetype='dotted', col='red') +
  geom_vline(xintercept = 2.6, linetype='dotted', col='red') +
  annotate("text", x = 2.6, y = 2.3, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(t3, (logval_bmi > 3.3 | logpval_crp > 3.3) & !grepl("^X", Variable)), 
                  aes(x=logval_bmi, y=logpval_crp, label=Variable), size = 3.5) +
  labs(color = "FDR significant")



plot_for_corrected <- ggplot(t3, aes(x=logval_bmi, y=logpval_crp)) +
  geom_point(aes(color = logval_bmi > 3.3 | logpval_crp > 3.3), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 6) +
  ylim(0, 6) +
  ggtitle('Differences in metabolomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 3.3, linetype='dotted', col='red') +
  geom_vline(xintercept = 3.6, linetype='dotted', col='red') +
  annotate("text", x = 3.3, y = 3.3, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(t3, (logval_bmi > 3.3 | logpval_crp > 3.3) & !grepl("^X", Variable)), 
                  aes(x=logval_bmi, y=logpval_crp, label=Variable), size = 3.5) +
  labs(color = "FDR significant")


#============================================================================
#TOTAL MODEL OF METABOLITES
clina <- clin_dat_at_diagnosis_clean %>% select(id, bs_bmi, cbt_inflammation_crp_mgpl)
clina <- na.omit(clina)

clinb <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex)
clinb <- na.omit(clinb)

follow_up_data <- readRDS("~/filepath/follow_up_data.rds")
fu <- follow_up_data %>% select(id, id_cohort) %>% unique()
clina <- merge(clina, fu, by='id')
clina <- na.omit(clina)
#now merge these 617 patients with metabolomics
clina$log_crp <- log(clina$cbt_inflammation_crp_mgpl)
clina <- merge(clina, clinb, by='id')
meta <- merge(clina, metablomics_normalisd_z_score_metabolomics, by.x='id_cohort', by.y='SUBJECT_ID')
#202 patients retained

numeric_cols <- meta %>% select(-id, -bs_bmi, -cbt_inflammation_crp_mgpl, -log_crp, - id_cohort, -age_diagnosis, -sex) %>% select(where(is.numeric))

library(broom)

results <- data.frame()

#get loop
for (col_name in colnames(numeric_cols)) {
  model <- lm(meta[[col_name]] ~ log_crp + bs_bmi + age_diagnosis + sex, data = meta)
  model_summary <- summary(model)$coefficients
  log_crp_coef <- model_summary["log_crp", "Estimate"]
  log_crp_pval <- model_summary["log_crp", "Pr(>|t|)"]
  bs_bmi_coef <- model_summary["bs_bmi", "Estimate"]
  bs_bmi_pval <- model_summary["bs_bmi", "Pr(>|t|)"]
  results <- rbind(results, data.frame(
    variable = col_name,
    log_crp_coef = log_crp_coef,
    log_crp_pval = log_crp_pval,
    bs_bmi_coef = bs_bmi_coef,
    bs_bmi_pval = bs_bmi_pval
  ))
}

results$fdr_crp <- p.adjust(results$log_crp_pval, method='fdr')
results$fdr_bmi <- p.adjust(results$bs_bmi_pval, method='fdr')

results$log10crp <- -log10(results$log_crp_pval)
results$log10bm <- -log10(results$bs_bmi_pval)

r1 <- results %>% filter(fdr_crp <0.05)
r2 <- results %>% filter(fdr_bmi <0.05)

getwd()
write.csv2(results, file='metabolomics_linear_model_adjusted.csv')

#now plot!
ggplot(results, aes(x=log10bm, y=log10crp)) +
  geom_point(aes(color = log10bm > 3 | log10crp > 4.3), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 6) +
  ylim(0, 6) +
  ggtitle('Differences in metabolomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 4.3, linetype='dotted', col='red') +
  geom_vline(xintercept = 3, linetype='dotted', col='red') +
  annotate("text", x = 3.3, y = 4.1, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(results, (log10bm > 3 | log10crp > 4.3) & !grepl("^X", variable)), 
                  aes(x=log10bm, y=log10crp, label=variable), size = 3.5) +
  labs(color = "FDR significant")

#=======================================================================================
#move on to proteomics
#we know these are normal
proteomics_correct <- readRDS("~/filepath/proteomics_correct.rds")
proteomics_correct <- proteomics_correct[,-c(1,3:8,10:64)]
prot <- proteomics_correct %>% filter(Event == 'PATIENT_1')
prot[,3:5286] <- lapply(prot[,3:5286], as.numeric)
prot <- as.data.frame(prot)

prot2 <- merge(weight, prot, by.x='id', by.y='OpenClinica.ID')

prot2$weight <- factor(prot2$weight, levels = c("underweight", "normal weight", "pre-obesity", 
                                                "obesity class-I", "obesity class-II", "obesity class-III"), 
                       ordered = TRUE)
#==============================================================
#quick look at FSTL3

prot2a <- prot2 %>% filter(!is.na(weight))

# Perform ANOVA
anova_result <- aov(FSTL3 ~ factor(weight), data = prot2a)
anova_pval <- summary(anova_result)[[1]]$`Pr(>F)`[1]  # Extract the p-value as a numeric value


# Create boxplot with ANOVA p-value annotation
ggplot(prot2a, aes(x = factor(weight), y = FSTL3, fill=weight)) +
  geom_boxplot() +
  labs(x = "Weight", y = "FSTL3 level") +
  theme_minimal() +
  ggtitle("Boxplot of FSTL3 by BMI group") +
  annotate("text", x = Inf, y = Inf, label = paste("ANOVA p =", formatC(anova_pval, format = "e", digits = 2)), 
           hjust = 1.1, vjust = 1.5, size = 4, color = "blue")


#===============================================================


numeric_cols <- prot2 %>% select(-id, -weight, -Event) %>% select(where(is.numeric))

# do ANOVA
anova_results <- lapply(names(numeric_cols), function(col_name) {
  formula <- as.formula(paste(col_name, "~ weight"))
  
  # Perform linear regression to obtain slope and intercept
  lm_result <- lm(formula, data = prot2)
  
  # Extract p-value from the ANOVA on the linear model
  p_value <- summary(aov(lm_result))[[1]][["Pr(>F)"]][1]
  
  # Extract intercept and slope
  coefficients <- coef(lm_result)
  intercept <- coefficients[1]
  slope <- coefficients[2]
  
  # Combine into a data frame
  data.frame(
    Variable = col_name,
    Intercept = intercept,
    Slope = slope,
    p_value = p_value
  )
})


# Combine results into a single dataframe
anova_summary_df <- do.call(rbind, anova_results)

# Adjust for FDR using p.adjust
anova_summary_df$p_value_fdr <- p.adjust(anova_summary_df$p_value, method = "fdr")

anova_summary_df$FDR <- p.adjust(anova_summary_df$p_value, method='fdr')
a2 <- anova_summary_df %>% filter(FDR <0.05)

Proteomics_names_conversion <- readRDS("~/filepath/Proteomics_names_conversion.rds")

prot_conv <- Proteomics_names_conversion
prot_conv$target_true <- gsub("[, /()\\-]", ".", prot_conv$Target)
#add numbering:
prot_conv$target_true <- ave(prot_conv$target_true, prot_conv$target_true, FUN = function(x) ifelse(seq_along(x) == 1, x, paste0(x, ".", seq_along(x) - 1)))
#select relevant data
pc <- prot_conv %>% select(EntrezGeneSymbol, target_true)
write_rds(pc, file='proteomics_easy_conversion_key.rds')


#order the factor
numeric_cols <- prot3 %>% select(-id, -above_5) %>% select(where(is.numeric))

# Ensure `above_5` is a factor
prot3$above_5 <- as.factor(prot3$above_5)

#do t-test
t_test_results <- lapply(names(numeric_cols), function(col_name) {
  # Check if both groups in 'above_5' have more than one unique value in the numeric column
  if(length(unique(prot3[[col_name]][prot3$above_5 == unique(prot3$above_5)[1]])) > 1 && 
     length(unique(prot3[[col_name]][prot3$above_5 == unique(prot3$above_5)[2]])) > 1) {
    
    # Create the formula for the t-test
    formula <- as.formula(paste(col_name, "~ above_5"))
    
    # Perform t-test
    t_test_result <- t.test(formula, data = prot3)
    
    # Extract p-value
    p_value <- t_test_result$p.value
    
    # Calculate the slope as the mean difference between the two groups
    group_means <- tapply(prot3[[col_name]], prot3$above_5, mean, na.rm = TRUE)
    slope <- diff(group_means)  # difference between the means of the two groups
    
    # Return results as a data frame
    return(data.frame(Variable = col_name, p_value = p_value, Slope = slope))
  } else {
    # If there's no variability, return NA for p-value and slope
    return(data.frame(Variable = col_name, p_value = NA, Slope = NA))
  }
})

# Combine results into a single dataframe
t_test_summary_df <- do.call(rbind, t_test_results)

# Adjust p-values using FDR (False Discovery Rate) correction
t_test_summary_df$FDR <- p.adjust(t_test_summary_df$p_value, method='fdr')

# Filter for significant results (FDR < 0.05)
significant_t_tests <- t_test_summary_df %>% filter(FDR < 0.05)

#test overlap
abc123 <- merge(a2, significant_t_tests, by='Variable')
#182 overlap ~50%
#=============================================================
#plot Chris' requested plot
#get plot for Chris
ta <- t_test_summary_df 
names(ta) <- c('Variable', 'CRP_pval', 'CRP_slope', 'CRP_FDR')

tb <- anova_summary_df
names(tb) <- c('Variable', 'bmi_intercept', 'BMI_slope', 'BMI_pval', 'fdr2','BMI_FDR')

tab <- merge(ta, tb, by='Variable')
tab$significance <- ifelse(tab$BMI_FDR < 0.05 & tab$CRP_FDR <0.05, 'both', NA)
tab$significance <- ifelse(tab$BMI_FDR < 0.05 & tab$CRP_FDR >=0.05, 'BMI only', tab$significance)
tab$significance <- ifelse(tab$BMI_FDR >=0.05 & tab$CRP_FDR <0.05, 'CRP only', tab$significance)
tab$significance <- ifelse(tab$BMI_FDR >=0.05 & tab$CRP_FDR >=0.05, 'none', tab$significance)

tab <- tab %>% filter(!is.na(significance))

ggplot(tab, aes(x=BMI_slope, y=CRP_slope, colour=significance)) +
  geom_point(size = 2, shape = 19) +  
  theme_bw() +
  xlim(-3.3, 3.3) +
  ylim(-3.3, 3.3) +
  ggtitle('Differences in proteomic signature') +
  xlab('Slope for BMI groups') +
  ylab('Slope for CRP groups') +
  labs(color = "FDR significance") +
  scale_color_manual(values = c("CRP only" = "red", "none" = "black", "BMI only" = "blue", "both" = "purple"))

#####################################################
unique(significant_t_tests$Variable)
#================================================================================
#make plot

#now plot
t1 <- t_test_summary_df 
t1$logpval_crp <- -log10(t1$p_value)

t2 <- anova_summary_df
t2$logval_bmi <- -log10(t2$p_value)

t3 <- merge(t1, t2, by='Variable')



#now plot
library(ggrepel)

t3 <- t3 %>% filter(!is.na(logpval_crp))

plot_scaled <- ggplot(t3, aes(x=logval_bmi, y=logpval_crp)) +
  geom_point(aes(color = logval_bmi > 2.6 | logpval_crp > 2.3), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 8) +
  ylim(0, 8) +
  ggtitle('Differences in proteomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 2.6, linetype='dotted', col='red') +
  geom_vline(xintercept = 2.39, linetype='dotted', col='red') +
  annotate("text", x = 2.39, y = 2.6, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(t3, (logval_bmi > 2.39 | logpval_crp > 2.6)), 
                  aes(x=logval_bmi, y=logpval_crp, label=Variable), size = 3.5) +
  labs(color = "FDR significant")


#now replicate this analysis for the genes of interest
#first select genes of interest:
library(readxl)
genes <- read_excel('/filepath/hpc-work/CRP analyses/Genes_of_interest_inflammation.xlsx')
interst <- c(genes$gene_name)

gene2 <- merge(genes, Proteomics_names_conversion, by.x='gene_name', by.y='EntrezGeneSymbol')
#91 out of the 94 genes have proteomics! 

gene2$target_new <- gsub("[- ]", ".", gene2$Target)
#NOTE TO SELF: ALSO REMOVE / AND , IF NEEDED LATER ON!
#now check which bit are not in gene2 
gene2_filtered <- gene2[!gene2$target_new %in% ab1$Variable, ]
gene2_filtered <- gene2_filtered %>% select(Target, target_new)
#4 need adjustment, some also overlap (as genes have 2 proteins!)

#change the following
gene2$target_new <- ifelse(gene2$target_new == '6Ckine', 'X6Ckine', gene2$target_new)
gene2$target_new <- ifelse(gene2$target_new == 'Fractalkine/CX3CL.1', 'Fractalkine.CX3CL.1', gene2$target_new)
gene2$target_new <- ifelse(gene2$target_new == 'CXCL16,.soluble', 'CXCL16..soluble', gene2$target_new)
gene2$target_new <- ifelse(gene2$target_new == 'Fas.ligand,.soluble', 'Fas.ligand..soluble', gene2$target_new)


ab1 <- anova_summary_df %>% filter(Variable %in% gene2$target_new)
ab2 <- t_test_summary_df %>% filter(Variable %in% gene2$target_new)
#so 84/94 genes have proteomics! 

#re-do FDR

ab1$fdr_bmi <- p.adjust(ab1$p_value, method='fdr')
ab1$pbmi <- -log10(ab1$p_value)
ab2$fdr_CRP <- p.adjust(ab2$p_value, method='fdr')
ab1$pcrp <- -log10(ab2$p_value)

ab3 <- merge(ab1, ab2, by='Variable')

#now plot
ggplot(ab3, aes(x=pbmi, y=pcrp)) +
  geom_point(aes(color = pbmi > 2.6 | pcrp > 2.3), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 11.5) +
  ylim(0, 11.5) +
  ggtitle('Differences in proteomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 2.2, linetype='dotted', col='red') +
  geom_vline(xintercept = 2.5, linetype='dotted', col='red') +
  annotate("text", x = 2.5, y = 2.2, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(ab3, (pbmi > 2.5 | pcrp > 2.2)), 
                  aes(x=pbmi, y=pcrp, label=Variable), size = 3.5) +
  labs(color = "FDR significant")



#==========================================================================================
#re-do joint model
#TOTAL MODEL OF METABOLITES
clina <- clin_dat_at_diagnosis_clean %>% select(id, bs_bmi, cbt_inflammation_crp_mgpl)
clina <- na.omit(clina)

clinb <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex)
clinb <- na.omit(clinb)

follow_up_data <- readRDS("~/projects/Data_for_cynapse/clinical_data/follow_up_data.rds")
fu <- follow_up_data %>% select(id, id_cohort) %>% unique()
clina <- merge(clina, fu, by='id')
clina <- na.omit(clina)
#now merge these 617 patients with proteomics
clina$log_crp <- log(clina$cbt_inflammation_crp_mgpl)
clina <- merge(clina, clinb, by='id')

p2 <- merge(clina, prot, by.x='id', by.y='OpenClinica.ID')
#319 patients retained

numeric_cols <- p2 %>% select(-id, -bs_bmi, -cbt_inflammation_crp_mgpl, -log_crp, - id_cohort, -age_diagnosis, -sex) %>% select(where(is.numeric))

library(broom)
results <- data.frame()

#get loop
for (col_name in colnames(numeric_cols)) {
  model <- lm(p2[[col_name]] ~ log_crp + bs_bmi+ age_diagnosis + sex, data = p2)
  model_summary <- summary(model)$coefficients
  log_crp_coef <- model_summary["log_crp", "Estimate"]
  log_crp_pval <- model_summary["log_crp", "Pr(>|t|)"]
  bs_bmi_coef <- model_summary["bs_bmi", "Estimate"]
  bs_bmi_pval <- model_summary["bs_bmi", "Pr(>|t|)"]
  results <- rbind(results, data.frame(
    variable = col_name,
    log_crp_coef = log_crp_coef,
    log_crp_pval = log_crp_pval,
    bs_bmi_coef = bs_bmi_coef,
    bs_bmi_pval = bs_bmi_pval
  ))
}

results$fdr_crp <- p.adjust(results$log_crp_pval, method='fdr')
results$fdr_bmi <- p.adjust(results$bs_bmi_pval, method='fdr')

results$log10crp <- -log10(results$log_crp_pval)
results$log10bm <- -log10(results$bs_bmi_pval)

write.csv2(results, file='proteomic_linear_models_bmi_crp.csv')

r1 <- results %>% filter(fdr_crp <0.05)
r1 <- merge(r1, pc, by.x='variable', by.y='target_true')
r2 <- results %>% filter(fdr_bmi <0.05)
r2 <- merge(r2, pc, by.x='variable', by.y='target_true')

#now plot!
ggplot(results, aes(x=log10bm, y=log10crp)) +
  geom_point(aes(color = log10bm > 2.6 | log10crp > 4), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 10.5) +
  ylim(0, 10.5) +
  ggtitle('Differences in proteomic signature') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 4, linetype='dotted', col='red') +
  geom_vline(xintercept = 2.6, linetype='dotted', col='red') +
  annotate("text", x = 3.3, y = 4.1, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(results, (log10bm > 2.6 | log10crp > 4)), 
                  aes(x=log10bm, y=log10crp, label=variable), size = 3.5) +
  labs(color = "FDR significant")


#=======================================================================================
#repeat proteomics in targeted approach
p2 <- merge(clina, prot, by.x='id', by.y='OpenClinica.ID')
#319 patients retained


numeric_cols <- p2 %>% select(-id, -bs_bmi, -cbt_inflammation_crp_mgpl, -log_crp, -id_cohort, -age_diagnosis, -sex) %>% select(where(is.numeric)) %>%  select(any_of(gene2$target_new))


library(broom)
results <- data.frame()

#get loop
for (col_name in colnames(numeric_cols)) {
  model <- lm(p2[[col_name]] ~ log_crp + bs_bmi+ age_diagnosis + sex, data = p2)
  model_summary <- summary(model)$coefficients
  log_crp_coef <- model_summary["log_crp", "Estimate"]
  log_crp_pval <- model_summary["log_crp", "Pr(>|t|)"]
  bs_bmi_coef <- model_summary["bs_bmi", "Estimate"]
  bs_bmi_pval <- model_summary["bs_bmi", "Pr(>|t|)"]
  results <- rbind(results, data.frame(
    variable = col_name,
    log_crp_coef = log_crp_coef,
    log_crp_pval = log_crp_pval,
    bs_bmi_coef = bs_bmi_coef,
    bs_bmi_pval = bs_bmi_pval
  ))
}

results$fdr_crp <- p.adjust(results$log_crp_pval, method='fdr')
results$fdr_bmi <- p.adjust(results$bs_bmi_pval, method='fdr')

results$log10crp <- -log10(results$log_crp_pval)
results$log10bm <- -log10(results$bs_bmi_pval)


#now plot!
ggplot(results, aes(x=log10bm, y=log10crp)) +
  geom_point(aes(color = log10bm > 2.7 | log10crp > 2.7), size = 2, shape = 19) +  
  scale_color_manual(values = c("black", "orange")) +  # Black for non-significant, orange for significant points
  theme_bw() +
  xlim(0, 10.5) +
  ylim(0, 10.5) +
  ggtitle('Differences in proteomic signature - age & sex adjusted') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 2.7, linetype='dotted', col='red') +
  geom_vline(xintercept = 2.7, linetype='dotted', col='red') +
  annotate("text", x = 6.5, y = 2.6, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(results, (log10bm > 2.7 | log10crp > 2.7)), 
                  aes(x=log10bm, y=log10crp, label=variable), size = 3.5) +
  labs(color = "FDR significant")


