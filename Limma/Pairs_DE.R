# library paks 
library(tidyverse)
library(limma)
library(GEOquery)

# load data 
GSE194314 <- getGEO('GSE194314', destdir=".",getGPL = F)
exprSet <- exprs(GSE194314[[1]]) 
boxplot(log2(exprSet))

# 标准化+log2
exprSet <- normalizeBetweenArrays(exprSet) %>% 
  log2(.)
boxplot(exprSet)

# 获取phenoData
pdata <- pData(GSE194314[[1]])

# 分组factor
individuals <- factor(unlist(lapply(pdata$characteristics_ch1.1,function(x) strsplit(as.character(x),":")[[1]][2])))

treatment <- unlist(lapply(pdata$characteristics_ch1.2,function(x) strsplit(as.character(x),":")[[1]][2]))

treatment <- factor(treatment,levels = unique(treatment))

# design matrix
# non-paired
design_non_paried <- model.matrix(~ 0 + treatment)
colnames(design_non_paried) <- c("Control","anti-BTLA")

fit1 <- lmFit(exprSet,design_non_paried)
fit1 <- eBayes(fit1)

# DE analysis
Diff_non_paired <- topTable(fit1,
                               adjust = 'BH',
                               coef =  "anti-BTLA",
                               n = Inf,
                               #p.value = 0.05
                               )

# paired
design_paried <- model.matrix(~ individuals + treatment)
summary(individuals)
summary(treatment)

fit2 <- lmFit(exprSet,design_paried)
fit2 <- eBayes(fit2)

# DE analysis
Diff_paired <- topTable(fit2,
                           adjust = 'BH',
                           coef = "treatment anti-BTLA",
                           n = Inf)

# plot
library(EnhancedVolcano)
library(patchwork)

p1 <- EnhancedVolcano(Diff_non_paired,
                      lab = rownames(Diff_non_paired),
                      x = 'logFC',
                      y =  'P.Value',
                      title = 'non_paired',
                      pointSize = 3.0,
                      labSize = 6.0,
                      legendPosition = 'right',
                      pCutoff = 0.05,
                      FCcutoff = 1)


p2 <- EnhancedVolcano(Diff_paired,
                      lab = rownames(Diff_paired),
                      x = 'logFC',
                      y =  'P.Value',
                      title = 'paired',
                      pointSize = 3.0,
                      labSize = 6.0,
                      legendPosition = 'right',
                      pCutoff = 0.05,
                      FCcutoff = 1)

p1 + p2
