# load packages
library(ggplot2)
library(e1071)
library(openintro)
library(dplyr)
library(nlme)
library(sqldf)
library(AICcmodavg)
library(grid)
library(gridExtra)
library(lattice)
library(plotrix)
library(gplots)
library(pheatmaps)
library(tidyr)
library(stringr)
library(ggdendro)
library(reshape2)
library(grid)
library(circlize)
library(olsrr)

# import data
data = read.csv(file="metab_cleaned.csv")

# removing the first column that R created when I saved it
data = data[,c(2:40)]

#############################################################################
# lm to determine relationship between follicle size and metabolite peak area
#############################################################################

# creating an empty dataframe to rbind my pvalues to:
lm_pvalues = data.frame(matrix(ncol = 2, nrow = 0))
col.names = c("metabolite", "pvalue")
colnames(lm_pvalues) = col.names

for(i in 2:ncol(data)) {      
    results = lm(data[,i]~data$g2_measure)
    s = summary(results)
    my.p = s$coefficients[2,4]
    new_row = data.frame(colnames(data[i]),my.p)
    lm_pvalues = rbind(lm_pvalues, new_row)
}
colnames(lm_pvalues) = col.names

# adjusting the pvalues using bonferroni and fdr
lm_pvalues$bonf = p.adjust(lm_pvalues$pvalue, method = "bonferroni")
lm_pvalues$fdrtest = p.adjust(lm_pvalues$pvalue, method = "fdr")

# outputing as a csv for all pvalues:
write.csv(lm_pvalues, file = "all_pvalues.csv")

# pulling the pvalues for the metabolites that are fdr significant
fdr_sig_pvalues = sqldf('select * from lm_pvalues where fdrtest <= 0.05')

# saving to a csv
write.csv(fdr_sig_pvalues, file = "fdr_sig_pvalues.csv")

# Only 18 significant metabolites
