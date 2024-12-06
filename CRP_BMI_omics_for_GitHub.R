#CRP and proteomics joint script
library(tidyverse)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(data.table)

#quick analysis
#EdB started here on 14-02-2024
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/filepath/PhD/Projects")
getwd()

RNA_expression_of_inflammatory_genes <- readRDS("C:/filepath/RNA_expression_of_inflammatory_genes.rds")
lapply(RNA_expression_of_inflammatory_genes[,2:86], hist)
RNA_expression_of_inflammatory_genes$above_5 <- as.factor(RNA_expression_of_inflammatory_genes$above_5)
proteomics_forSave_crp_bmi <- readRDS("C:/filepath/proteomics_forSave_crp_bmi.rds")

#now log-transform the data
rna1 <- RNA_expression_of_inflammatory_genes
rna2 <- rna1

rna1[,2:86] <- rna1[,2:86] +10
rna1[,2:86] <- lapply(rna1[,2:86], log)
lapply(rna1[,2:86], hist)
#data is much more normal but still some clearly bimodal distributions present....
rna <- rna1 %>% pivot_longer(2:86, names_to = 'gene', values_to = 'expression_level')
df <- as.data.table(rna)
res1 <- df[, t.test(data=.SD, expression_level ~ above_5), by=gene]
res1 <- res1 %>% select(gene, p.value)
res1 <- unique(res1)
res1$qval <- p.adjust(res1$p.value, method = 'fdr')


#also do wilcox test
rna1$above_5 <- as.factor(rna1$above_5)

mann_whit = function(myvar) {
  a <- wilcox.test(rna1[[myvar]] ~ rna1[['above_5']], data=rna1)
  print(a$p.value)
}

#check function
mann_whit('BDNF')

#function works, now into lapply
vars <- c(colnames(rna1[,2:86]))

pval <- lapply(vars, mann_whit)
pval2 <- as.data.frame(pval)
pval2 <- t(pval2)


#############################################################################
mann_whit <- function(myvar) {
  test_result <- wilcox.test(rna1[[myvar]] ~ rna1[['above_5']], data = rna1)
  result <- data.frame(
    variable = myvar,
    p_value = test_result$p.value)
  
  return(result)
}

results <- do.call(rbind, lapply(vars, mann_whit))
print(results)

#=============================================================================
#compare for weight groups
inflammatory_genes_and_adipokines <- readRDS("C:/filepath/inflammatory_genes_and_adipokines.rds")
#now log-transform the data
rna1 <- inflammatory_genes_and_adipokines

rna1[,2:101] <- rna1[,2:101] +10
rna1[,2:101] <- lapply(rna1[,2:101], log)
lapply(rna1[,2:101], hist)
#data is much more normal but still some clearly bimodal distributions present....

rna <- rna1 %>% pivot_longer(2:101, names_to = 'gene', values_to = 'expression_level')

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/filepath/v4_clean_clinical_data_first_visit_15Dec23.rds")
bmi_df <- v4_clean_clinical_data_first_visit_15Dec23
bmi_df$weight <- ifelse(bmi_df$bs_bmi <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=18.5 & bmi_df$bs_bmi <25, 'normal weight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=25 & bmi_df$bs_bmi <30, 'pre-obesity', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=30 &  bmi_df$bs_bmi <35, 'obesity class-I', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=35 &  bmi_df$bs_bmi <40, 'obesity class-II', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=40, 'obesity class-III', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)


bmi_df <- bmi_df %>% select(id, weight)
bmi_df <- bmi_df %>% filter(!is.na(weight))
rna_bmi <- merge(rna, bmi_df, by.x='OpenClinicaID', by.y='id')
#n=356 for BMI groups

#repeat with kruskal
krusk_all_genes = function(myvar) {
  bmi2 <- rna_bmi %>% filter(gene == myvar) 
  krusk <- kruskal.test(expression_level ~ weight, data=bmi2)
  return(tidy(krusk))
}

#test
krusk_all_genes('ANGPT2')
#this works

#now get a function with lapply
myvar = unique(rna_bmi$gene)
output = lapply(myvar, krusk_all_genes)
names(output) = myvar

combined_output = dplyr::bind_rows(output, .id = "gene")

bmi_res <- combined_output[,c(1,3)]
crp_res <- results

names(bmi_res) <- c('gene', 'pval_bmi')
names(crp_res) <- c('gene', 'pval_crp')

#==================================================================
genes <- merge(bmi_res, crp_res, by='gene')

#now scale the pvalues and plot
genes$FDR_adjusted_for_CRP <- p.adjust(genes$pval_crp, method='fdr')
genes$log10CRP <-  -log10(genes$pval_crp)
genes$FDR_adjusted_for_BMI <- p.adjust(genes$pval_bmi, method='fdr')
genes$log10BMI <-  -log10(genes$pval_bmi)

#now plot
plot_scaled <- ggplot(genes, aes(x=log10BMI, y=log10CRP)) +
  geom_point(aes(color = ifelse(log10CRP > 2 | log10BMI > 2, "orange", "black"))) +
  scale_color_identity() + 
  theme_bw() +
  xlim(0, 4.75) +
  ylim(0,4.75) +
  ggtitle('Differences in inflammatory gene expression') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 2, linetype='dotted', col='red') +
  geom_vline(xintercept = 2, linetype='dotted', col='red') +
  annotate("text", x = 3.5, y = 2, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(genes, log10CRP >2 | log10BMI >2), aes(x=log10BMI, y=log10CRP, label=gene))

#=====================================================================================
#now repeat this for the proteomics
crp_high_low_first_or_diagnostic_visit_v131Jan24 <- readRDS("C:/filepath/crp_high_low_first_or_diagnostic_visit_v131Jan24 (1).rds")

crp <- crp_high_low_first_or_diagnostic_visit_v131Jan24[,c(1,7)]
crp$above_5 <- as.factor(crp$above_5)

prot <- proteomics_forSave_crp_bmi

p1 <- merge(crp, prot, by.x='id', by.y='OpenClinica.ID')
#387 patients
unique(p1$id)
p2 <- merge(bmi_df, prot, by.x='id', by.y='OpenClinica.ID')
#441 patients 
unique(p2$id)

#===================================================
#now test CRP first
mann_whit <- function(myvar) {
  test_result <- wilcox.test(p1[[myvar]] ~ p1[['above_5']], data = p1)
  result <- data.frame(
    variable = myvar,
    p_value = test_result$p.value)
  
  return(result)
}

vars <- c(colnames(p1[,3:86]))
results <- do.call(rbind, lapply(vars, mann_whit))
print(results)

#also do kruskal for BMI
pbmi <- p2 %>% pivot_longer(3:86, names_to = 'protein', values_to = 'value')

#repeat with kruskal
krusk_all_genes = function(myvar) {
  bmi2 <- pbmi %>% filter(protein == myvar) 
  krusk <- kruskal.test(value ~ weight, data=bmi2)
  return(tidy(krusk))
}

#now get a function with lapply
myvar = unique(pbmi$protein)
output = lapply(myvar, krusk_all_genes)
names(output) = myvar

combined_output = dplyr::bind_rows(output, .id = "protein")

bmi_res <- combined_output[,c(1,3)]
crp_res <- results

names(bmi_res) <- c('protein', 'pval_bmi')
names(crp_res) <- c('protein', 'pval_crp')

#==================================================================
proteins <- merge(bmi_res, crp_res, by='protein')

#now scale the pvalues and plot
proteins$FDR_adjusted_for_CRP <- p.adjust(proteins$pval_crp, method='fdr')
proteins$log10CRP <-  -log10(proteins$pval_crp)
proteins$FDR_adjusted_for_BMI <- p.adjust(proteins$pval_bmi, method='fdr')
proteins$log10BMI <-  -log10(proteins$pval_bmi)

#now plot
plot_protein <- ggplot(proteins, aes(x=log10BMI, y=log10CRP)) +
  geom_point(aes(color = ifelse(log10CRP > 2.4 | log10BMI > 2, "orange", "black"))) +
  scale_color_identity() + 
  theme_bw() +
  xlim(0, 20) +
  ylim(0,15) +
  ggtitle('Differences in inflammatory protein level') +
  xlab('-log10 P-values between BMI groups') +
  ylab('-log10 P-values between CRP groups') +
  geom_abline() +
  geom_hline(yintercept = 2.4, linetype='dotted', col='red') +
  geom_vline(xintercept = 2, linetype='dotted', col='red') +
  annotate("text", x = 10, y = 2.5, label = "FDR significance threshold", vjust = -0.5, size=3) +
  theme(aspect.ratio = 1) +
  geom_text_repel(data=subset(proteins, log10CRP >2.4 | log10BMI >2), aes(x=log10BMI, y=log10CRP, label=protein))

Proteomics_names_conversion <- readRDS("~/filepath/Omics/Proteomics_names_conversion.rds")
proteomics_easy_conversion_key <- readRDS("~/filepath/Omics/proteomics_easy_conversion_key.rds")

protnames <- merge(Proteomics_names_conversion[,c(4,8)], proteomics_easy_conversion_key, by='EntrezGeneSymbol')
protnames <- protnames %>% filter(EntrezGeneSymbol != '')

g1 <- genes[,c(1,4,6)]
p1 <- proteins[,c(1,4,6)]
p1 <- merge(p1, protnames, by.x='protein', by.y='target_true')
p1 <- p1[,-1]
p1<- unique(p1)

a1 <- merge(g1, p1, by.x='gene', by.y='EntrezGeneSymbol', all=T)

names(a1)<- c('gene', 'FDR adjusted P-value for CRP groups in transcriptomics', 'FDR adjusted P-value for BMI groups in transcriptomics', 'FDR adjusted P-value for CRP groups in proteomics', 'FDR adjusted P-value for BMI groups in proteomics', 'protein')
a1 <- a1[,c(1,6,2:5)]
write.csv2(a1, file='transcriptomics_proteomic_differences.csv')
getwd()