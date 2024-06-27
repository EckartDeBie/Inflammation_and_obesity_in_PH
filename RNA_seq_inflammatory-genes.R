#quick analysis
#EdB started here on 14-02-2024
data_dir="~/PhD/Projects" 
getwd()
setwd("C:/Users/Gebruiker/Documents/PhD/Projects")
getwd()

RNA_expression_of_inflammatory_genes <- readRDS("C:/filepath.rds")
lapply(RNA_expression_of_inflammatory_genes[,2:86], hist)
RNA_expression_of_inflammatory_genes$above_5 <- as.factor(RNA_expression_of_inflammatory_genes$above_5)

#now log-transform the data
rna1 <- RNA_expression_of_inflammatory_genes
rna2 <- rna1

rna1[,2:86] <- rna1[,2:86] +10
rna1[,2:86] <- lapply(rna1[,2:86], log)
lapply(rna1[,2:86], hist)
#data is much more normal but still some clearly bimodal distributions present....
library(tidyverse)

rna <- rna1 %>% pivot_longer(2:86, names_to = 'gene', values_to = 'expression_level')
library(data.table)
df <- as.data.table(rna)
res1 <- df[, t.test(data=.SD, expression_level ~ above_5), by=gene]
res1 <- res1 %>% select(gene, p.value)
res1 <- unique(res1)
res1$qval <- p.adjust(res1$p.value, method = 'fdr')

#also try to generate a plot
library(rstatix)
library(ggpubr)
ggplot(rna, aes(x=above_5, y=expression_level, fill=above_5)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = 'free') +
  stat_compare_means(method = 't.test')

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

res1$non_parametric_pval <- pval2
res1$qnonp <- p.adjust(res1$non_parametric_pval, method='fdr')
pvals_fdr_crp_groups <- res1 %>% select(gene, non_parametric_pval)
names(pvals_fdr_crp_groups) <- c('gene', 'pval_for_CRP')

#get logFC ready
#1) get mean expression data
rna_nolog <- rna2 %>% pivot_longer(2:86, names_to = 'gene', values_to = 'expression_level')

mean_express <- rna_nolog %>% group_by(above_5, gene) %>% summarise(mean_expression = mean(expression_level, na.rm=T))
expression_spread <- mean_express %>% pivot_wider(names_from = above_5, values_from = mean_expression)
expression_spread$logFC <- expression_spread$yes/expression_spread$no
expression_spread$logFC2 <- log2(expression_spread$logFC)

#merge this with qvals
volc_pl <- merge(res1, expression_spread, by='gene')
volc_pl <- as.data.frame(volc_pl)
rownames(volc_pl) <- volc_pl$gene

library(EnhancedVolcano)
EnhancedVolcano(volc_pl,
                lab = rownames(volc_pl),
                x = 'logFC2',
                y = 'qval',
                xlim = c(-1, 1),
                title = 'volcano plot of high versus low CRP inflammatory genes',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

#select the significantly different ones
res_plot <- res1 %>% filter(qnonp <0.05)

res_plot <- rna %>% filter(gene %in% res_plot$gene)


ggplot(res_plot, aes(x=above_5, y=expression_level, fill=above_5)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = 'free') +
  stat_compare_means(method = "t.test") +
  scale_y_continuous(expand=expansion(mult = c(0, .1))) +
  ggtitle('Differential, PC corrected, TPM, of inflammatory genes between high/low CRP')

#also with adjusted colours
ggplot(res_plot, aes(x=above_5, y=expression_level, fill=above_5)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = 'free') +
  scale_y_continuous(expand=expansion(mult = c(0, .1))) +
  scale_fill_manual(values=c('darkblue', 'darkred')) +
  ggtitle('Differential, PC corrected, TPM, of inflammatory genes between high/low CRP')

library(pheatmap)
library(RColorBrewer)
col1 <- c('darkblue', 'darkred')

hm <- res_plot %>% select(above_5, gene, expression_level, SampleID)

hm2 <- hm[,-1]
hm2 <- hm2 %>% pivot_wider(names_from = SampleID, values_from = expression_level)

hm2 <- as.data.frame(hm2)
rownames(hm2) <- hm2$gene
hm2 <- hm2[-1]
hm2 <- as.matrix(hm2)

meta_dat <- hm[,c(1,4)]

#get annotation
levels(meta_dat$above_5)
meta_dat <- unique(meta_dat)
rownames(meta_dat) <- meta_dat$SampleID
meta_dat2 <- meta_dat %>% select(above_5)
meta_dat2 <- as.data.frame(meta_dat2)

pheatmap(hm2, annotation_col = meta_dat, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = 204)

#this is all very difficult --> let's change this
#calculate relative difference
hm2 <- hm[,1:4]

hm2 <- hm2 %>% pivot_wider(names_from = gene, values_from = expression_level)
hm2 <- hm2 %>% select(! SampleID)

hm2[,2:8] <- lapply(hm2[,2:8], log)
hm2[,2:8] <- lapply(hm2[,2:8], scale)

hm3 <- hm2 %>% pivot_longer(2:8, names_to = 'gene', values_to = 'scaled_expression_level')
hm4 <- hm3 %>% group_by(above_5, gene) %>% summarise(median_expression = median(scaled_expression_level))
hm4 <- hm4 %>% pivot_wider(values_from = median_expression, names_from = above_5)
hm4 <- as.data.frame(hm4)
rownames(hm4) <- hm4$gene
hm4 <- hm4 %>% select(! gene)

colnames(hm4)<- c('CRP <=5', 'CRP >5')
hm5 <- as.matrix(hm4)
pheatmap(hm5, cluster_rows = FALSE, cluster_cols = FALSE, legend=TRUE)


#=============================================================================
#compare for weight groups
#CRP modelling without clustering
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
library(broom)


inflammatory_genes_and_adipokines <- readRDS("C:/filepaths.rds")
lapply(inflammatory_genes_and_adipokines[,2:101], hist)

#now log-transform the data
rna1 <- inflammatory_genes_and_adipokines

rna1[,2:101] <- rna1[,2:101] +10
rna1[,2:101] <- lapply(rna1[,2:101], log)
lapply(rna1[,2:101], hist)
#data is much more normal but still some clearly bimodal distributions present....

rna <- rna1 %>% pivot_longer(2:101, names_to = 'gene', values_to = 'expression_level')

v4_clean_clinical_data_first_visit_15Dec23 <- readRDS("C:/filepath.rds")
bmi_df <- v4_clean_clinical_data_first_visit_15Dec23
bmi_df$weight <- ifelse(bmi_df$bs_bmi <18.5, 'underweight', NA)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=18.5 & bmi_df$bs_bmi <25, 'normal weight', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=25 & bmi_df$bs_bmi <30, 'pre-obesity', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=30 &  bmi_df$bs_bmi <35, 'obesity class-I', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=35 &  bmi_df$bs_bmi <40, 'obesity class-II', bmi_df$weight)
bmi_df$weight <- ifelse(bmi_df$bs_bmi >=40, 'obesity class-III', bmi_df$weight)
bmi_df$weight <- as.factor(bmi_df$weight)


bmi_df <- bmi_df %>% select(id, weight)

rna_bmi <- merge(rna, bmi_df, by.x='OpenClinicaID', by.y='id')

anova_all_genes = function(myvar) {
  bmi2 <- rna_bmi %>% filter(gene == myvar) 
  anov <- aov(expression_level ~ weight, data=bmi2)
  return(tidy(anov))
}

#test
anova_all_genes('ANGPT2')
#this works

#now get a function with lapply
myvar = unique(rna_bmi$gene)
output = lapply(myvar, anova_all_genes)
names(output) = myvar

combined_output = dplyr::bind_rows(output, .id = "gene")
#remove the residuals

comb_output_anova <- combined_output %>% filter(term == 'weight')
comb_output_anova$fdr <- p.adjust(comb_output_anova$p.value, method = 'fdr')
comb1 <- comb_output_anova %>% filter(fdr <0.05)
#nothing passes FDR threshold!

#repeat with kruskal

krusk_all_genes = function(myvar) {
  bmi2 <- rna_bmi %>% filter(gene == myvar) 
  anov <- kruskal.test(expression_level ~ weight, data=bmi2)
  return(tidy(anov))
}

#test
krusk_all_genes('ANGPT2')
#this works

#now get a function with lapply
myvar = unique(rna_bmi$gene)
output = lapply(myvar, krusk_all_genes)
names(output) = myvar

combined_output = dplyr::bind_rows(output, .id = "gene")
#remove the residuals

comb_output_krusk <- combined_output
comb_output_krusk$fdr <- p.adjust(comb_output_krusk$p.value, method = 'fdr')
comb1 <- comb_output_anova %>% filter(fdr <0.05)

bmi_genes_diff <- comb_output_krusk %>% select(gene, p.value)
names(bmi_genes_diff) <- c('gene', 'pval_for_BMI')

genes <- merge(pvals_fdr_crp_groups, bmi_genes_diff, by='gene')

#now scale the pvalues and plot
genes$FDR_adjusted_for_CRP <- -log10(genes$pval_for_CRP)
genes$FDR_adjusted_for_BMI <- -log10(genes$pval_for_BMI)

lapply(genes[,2:3], hist)

genes <- as.data.frame(genes)

#now plot
library(ggrepel)

plot_scaled <- ggplot(genes, aes(x=FDR_adjusted_for_BMI, y=FDR_adjusted_for_CRP)) +
  geom_point() +
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
  geom_text_repel(data=subset(genes, FDR_adjusted_for_CRP >2 | FDR_adjusted_for_BMI >2), aes(x=FDR_adjusted_for_BMI, y=FDR_adjusted_for_CRP, label=gene))

            