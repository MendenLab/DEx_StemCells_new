################################################################################
###### RNA-seq analysis with DESeq2 of Ferroptosis & Antioxidant effect ########
################################################################################

## libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
BiocManager::install("DESeq2")
install.packages("BiocManager")
install.packages("ggpubr")
library(ggplot2)
library(magrittr)
library(ggpubr)
library(readxl)
library(readr)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(plyr)
library(ggplot2)
library(ggrepel)

## load data
countdata <- rall_counts <- read_delim("all.counts.tsv",
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)
View(countdata)

gene_names <- countdata[,1] %>% pull()
countdata <- countdata[,-1]
rownames(countdata) <- gene_names

## preparation for differential expression
# definition of conditions
coldata = data.frame(rownames = names(countdata),
                     condition = rep(c("Day 0",
                                       "Day 10",
                                       "Day 20 B27+AO",
                                       "Day 20 B27-AO",
                                       "Day 20 B27-AO+Fer",
                                       "Day 20 B27-AO+DMSO"),
                                     each=3)) 

# change "+" and "-" to avoid issues
coldata$condition <- gsub( "[+]", "plus", coldata$condition )
coldata$condition <- gsub( "[-]", "minus", coldata$condition )

# creation of DESeq object
dds <- DESeqDataSetFromMatrix (countData = countdata, 
                               colData = coldata, 
                               design = ~ condition)

## Pre-processing - removing low counts
# 28 516 -> 20 150
keep <- rowSums(counts(dds)) >= ncol(dds)
dds <- dds[keep,]

## DESeq2 object 
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds) 
coldata$condition # check names for conditions

## Results per condition
res0_10 <- results(dds, contrast = c("condition",'Day 0','Day 10')) # 0 vs 10
res10_20plusAO <- results(dds, contrast = c("condition","Day 10","Day 20 B27plusAO"))  #10 vs 20+AO
res0_20plusAO <- results(dds, contrast = c("condition",'Day 0','Day 20 B27plusAO')) #0 vs 20+AO

res_20minusAO_20plusA0 <- results(dds, contrast = c("condition",'Day 20 B27minusAO','Day 20 B27plusAO'))  # 20-AO vs 20+AO
res_20minusAO_20minusA0plusDMSO <- results(dds, contrast = c("condition",'Day 20 B27minusAO','Day 20 B27minusAOplusDMSO'))   # 20-AO vs 20+AO
res_20minusAO_20minusA0plusFer <- results(dds, contrast = c("condition",'Day 20 B27minusAO','Day 20 B27minusAOplusFer'))  # 20-AO vs 20+AO
res_20plusAO_20minusA0plusDMSO <- results(dds, contrast = c("condition",'Day 20 B27plusAO','Day 20 B27minusAOplusDMSO'))  # 20-AO vs 20+AO
res_20plusAO_20minusA0plusFer <- results(dds, contrast = c("condition",'Day 20 B27plusAO','Day 20 B27minusAOplusFer'))   # 20-AO vs 20+AO
res_20minusA0plusDMSO_20minusA0plusFer <- results(dds, contrast = c("condition",'Day 20 B27minusAOplusDMSO','Day 20 B27minusAOplusFer'))  # 20-AO vs 20+AO

res10_20minusAO <- results(dds, contrast = c("condition","Day 10","Day 20 B27minusAO"))   # 20-AO vs 20+AO
res10_20minusA0plusDMSO <- results(dds, contrast = c("condition","Day 10","Day 20 B27minusAOplusDMSO"))  # 20-AO vs 20+AO
res10_20minusA0plusFer <- results(dds, contrast = c("condition","Day 10","Day 20 B27minusAOplusFer"))   # 20-AO vs 20+AO

## volcano plot 
# auxiliar function to modify results matrices for volcano plot
results_matrix_plot <- function(res_matrix){
  res_matrix <- as.data.frame(res_matrix)
  res_matrix <- res_matrix %>%
    mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                 log2FoldChange <= (-1) & padj <= 0.05 ~ "down",
                                 TRUE ~ "ns"),
           hgnc_symbol = rownames(.)) 
}

res0_10 <- results_matrix_plot(res0_10)
res10_20plusAO <- results_matrix_plot(res10_20plusAO)
res0_20plusAO <- results_matrix_plot(res0_20plusAO)
res_20minusAO_20plusA0 <- results_matrix_plot(res_20minusAO_20plusA0)
res_20minusAO_20minusA0plusDMSO <- results_matrix_plot(res_20minusAO_20minusA0plusDMSO)
res_20minusAO_20minusA0plusFer <- results_matrix_plot(res_20minusAO_20minusA0plusFer) 
res_20plusAO_20minusA0plusDMSO <- results_matrix_plot(res_20plusAO_20minusA0plusDMSO) 
res_20plusAO_20minusA0plusFer <- results_matrix_plot(res_20plusAO_20minusA0plusFer)
res_20minusA0plusDMSO_20minusA0plusFer <- results_matrix_plot(res_20minusA0plusDMSO_20minusA0plusFer)
res10_20minusAO <- results_matrix_plot(res10_20minusAO)
res10_20minusA0plusDMSO <- results_matrix_plot(res10_20minusA0plusDMSO)
res10_20minusA0plusFer <- results_matrix_plot(res10_20minusA0plusFer)

genes_interest <- c("RARB","BCO2")

# Volcano plot 
volcano_plot <- function(data_res,title_comp){
  results_gene_interes <- data_res %>%
    filter(hgnc_symbol %in% genes_interest)
  cols <- c("up" = "#dd0d0b", "down" = "#140101", "ns" = "grey") 
  sizes <- c("up" = 4, "down" = 4, "ns" = 2) 
  alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)
  ggplot(data = data_res,
         aes(x = log2FoldChange,
             y = -log10(padj))) + 
    geom_point(aes(colour = gene_type), 
               alpha = 0.7, 
               shape = 16,
               size = 3) + 
    geom_point(data = results_gene_interes,
               shape = 21,
               size = 4,aes(fill = gene_type)) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "grey48") + 
    annotate(geom="text", x=3.8, y=(-log10(0.05)) + 2, 
             label="FDR = 5%",color = "grey48") +
    geom_vline(xintercept = c(- 1, 1),
               linetype = "dashed",
               color = "grey48") +
    geom_text_repel(data = subset(data_res, (gene_type != "ns" & padj < 0.0001)),
                     aes(label = hgnc_symbol),
                    size = 4, # font size in the text labels
                    point.padding = 0, # additional padding around each point
                    min.segment.length = 0.08, # draw all line segments
                    max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                    box.padding = 0.4,
                    nudge_x = .6,
                    nudge_y = .4,
                    hjust = 0.5) + # additional padding around each text label) +
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) +
    scale_x_continuous(limits = c(-4, 4)) +
    labs(title = "",
         x = expression("log"[2]*"(fold-change)"),
         y = expression("-log"[10]*"(adjusted p-value)"),
         colour = "Differential \nExpression") +
    theme_classic() + # Select theme with a white background  
    guides(fill=FALSE) +
    theme(axis.text=element_text(size=14),
           axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          legend.text=element_text(size=14),
          legend.position = c(0.1, 0.9)) 
  
}

# volcano plot of B27-AO vs B27+AO
volcano_plot(res_20minusAO_20plusA0,"B27-AO vs B27+AO")

## Boxplots of RARB and BCO2
tcounts <- t(log2((counts(dds[genes_interest, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(genes_interest)+1):ncol(.))

# only +AO and -AO on day 20
tcounts_filter <- tcounts %>%
  filter(condition %in% c("Day 20 B27plusAO","Day 20 B27minusAO"))
  
tcounts_filter$condition <- factor(tcounts_filter$condition,
                           labels = c("B27-AO","B27+AO"))
tcounts_filter$condition <- factor(tcounts_filter$condition,levels = c("B27+AO","B27-AO"))

for (i in genes_interest) {
  p <- ggplot(filter(tcounts_filter, gene==i), aes(condition, expression, color=condition)) + 
    geom_boxplot() + 
    geom_point(position=position_jitterdodge(),
               alpha = 0.7, 
               size = 3) +
    scale_color_manual(values=c("B27+AO" = "#140101", "B27-AO" = "#dd0d0b")) +
    labs(x = "",
         y = "Expression (log normalized counts)",
         fill = "",
         title = i) +
    stat_compare_means(method = "t.test", aes(label = sprintf("p-value = %2.1e", as.numeric(..p.format..))),
                       label.x = 1.25, 
                       label.y = 8,
                       size = 5)+
    theme_classic2() + # Select theme with a white background  
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          plot.title = element_text(size = 13, face="bold",hjust = 0.5),
          legend.text=element_text(size=12),
          legend.position="none") 
  print(p)
}

## Time-plot of RARB and BCO2
tcounts_mean <- tcounts %>%
  dplyr::group_by(condition,gene) %>%
  dplyr::mutate(expression_mean = mean(expression)) %>%
  dplyr::select(condition,expression_mean,gene) %>%
  filter(condition %in% c("Day 0", "Day 10","Day 20 B27plusAO", "Day 20 B27minusAO")) %>%
  distinct() 

tcounts_mean <- rbind(data.frame(tcounts_mean %>% filter(condition != "Day 20 B27plusAO") %>% mutate(type_condition = "B27 -AO")),
                      data.frame(tcounts_mean %>% filter(condition != "Day 20 B27minusAO") %>% mutate(type_condition = "B27 +AO")))

tcounts_mean$condition <- factor(tcounts_mean$condition,
                                   labels = c("Day 0","Day 10","Day 20","Day 20"))
tcounts_mean$condition <- factor(tcounts_mean$condition,levels = c("Day 0","Day 10","Day 20"))

for (i in genes_interest) {
  p <- ggplot(filter(tcounts_mean, gene==i), aes(condition, expression_mean,group = type_condition)) +  
    geom_point(aes(color=type_condition),
               alpha = 0.7, 
               size = 3) +
    geom_line(aes(color=type_condition),
              alpha = 0.8) +
    scale_color_manual(values=c("B27 +AO" = "#140101", "B27 -AO" = "#dd0d0b")) +
    labs(x = "",
         y = "Mean Expression (log normalized counts)",
         color = "",
         title = i) +
    theme_minimal() + # Select theme with a white background  
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          plot.title = element_text(size = 13, face="bold",hjust = 0.5),
          legend.text=element_text(size=12)) 
  print(p)
}

tcounts_mean$type_condition <- ifelse(tcounts_mean$condition == "Day 0","Day 0",
                                      ifelse(tcounts_mean$condition == "Day 10","Day 10",tcounts_mean$type_condition))

for (i in genes_interest) {
  p <- ggplot(filter(tcounts_mean, gene==i), aes(x = condition,y= expression_mean,fill = type_condition)) +  
    # geom_point(aes(color=type_condition),
    #            alpha = 0.7, 
    #            size = 3) +
    geom_bar(stat = "identity",position = "dodge",
             alpha = 0.6) +
    scale_fill_manual(values=c("B27 +AO" = "#140101", "B27 -AO" = "#dd0d0b", "Day 0" = "grey48", "Day 10" = "grey48")) +
    labs(x = "",
         y = "Mean Expression (log normalized counts)",
         color = "",
         title = i) +
    theme_minimal() + # Select theme with a white background  
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          plot.title = element_text(size = 13, face="bold",hjust = 0.5),
          legend.text=element_text(size=12)) 
  print(p)
}