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

##################################################################
# ANOVA and lm to determine effects of cow traits on follicle size
##################################################################

# import data
data = read.csv(file = "metabolomics_2020.csv")
cow_info = read.csv(file = "cow_info.csv")
foll_info = read.csv(file = "foll_class_info.csv")

# pulling out only the 2020 cows
cow_2020 = sqldf('select * from cow_info where year = 2020')
foll_2020 =sqldf('select * from foll_info where year = 2020')

# pulling out only the cows that were used in the metabolome study
met_cows = sqldf('select * from cow_2020 where cow_id in (select cow from data)') # used a subquery within sqldf to generate a list of cow ids from the metabolome data and then selected cow info from 2020 based on that list)
met_foll = sqldf('select * from foll_2020 where cow_id in (select cow from data)') # repeated with foll size info

# creating a table with desired info from met_cows and met_foll:
met_data = sqldf('select c.cow_id as id, f.g2_measure, f.foll_class, c.age, c.bcs, c.wt2, c.days_pp_fa, c.g1_to_pg, c.pg_to_g2, c.g2_to_fa from met_cows as c left join met_foll as f where c.cow_id = f.cow_id')

# converting weight in pounds to weight in kg
met_data$wt2 = met_data$wt2/2.20462
met_data$pg_to_g2 = met_data$pg_to_g2*24


# removing outlier cows (unknown follicle size/aborted)
met_data = sqldf('select * from met_data where id is not 1767')
met_data = sqldf('select * from met_data where id is not 2055')
met_data = sqldf('select * from met_data where id is not 2778')

met_data = droplevels(met_data)

# writing to a csv so that I can skip all the above next time
write.csv(met_data, file = "met_data_cleaned.csv")

# for loop to do lm on the different cow traits relative to follicle size and add the p value to a dataframe
# creating an empty dataframe to append the pvalues to based on size only
pvalues = data.frame(matrix(ncol = 2, nrow = 0))
col.names = c("cow trait", "pvalue")
colnames(pvalues) = col.names
# running the loop
for(i in 4:5) {     
    results = aov(met_data[,i]~met_data$g2_measure)
    s = summary(results)
    my.p = summary(results)[[1]][1,5]
    new_row = data.frame(colnames(met_data[i]),my.p)
    pvalues = rbind(pvalues, new_row)
}
for(i in 6:ncol(met_data)) {
    results = lm(met_data[,i]~met_data$g2_measure)
    s = summary(results)
    my.p = s$coefficients[2,4]
    new_row = data.frame(colnames(met_data[i]),my.p)
    pvalues = rbind(pvalues, new_row)
}
# changing the column names
colnames(pvalues) = col.names

# exporting it into an excel file for future reference:
write.csv(pvalues, file = "cow_trait_pvalues.csv")


# creating plots for each lm for each of the cow traits vs follicle size
# age
met_data$age = as.factor(met_data$age)
my.p = round(pvalues$pvalue[1], digits = 3)
text_p=paste(toString("p="), toString(my.p), sep="")
jpeg("age.jpeg")
    age_plot = ggplot(met_data, aes(x=age, y=g2_measure)) + 
    geom_violin() + 
    geom_point(size = 3) + 
    xlab("Age (yrs)") + 
    ylab("Follicle Size (mm)") + 
    annotate(geom = 'text', label = text_p, x = 0.75, y = 17.8, hjust = 0, vjust = 1, cex = 8) + 
    ggtitle("A") + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 23, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05),
    axis.title = element_text(size = 23, color = "black", face = "bold")) + 
    stat_summary(fun.y=mean, geom="crossbar",width = 0.2, color = "#0072B2")
    print(age_plot)
dev.off()



# BCS
met_data$bcs = as.factor(met_data$bcs)
my.p = round(pvalues$pvalue[2], digits = 3)
text_p=paste(toString("p="), toString(my.p), sep="")
jpeg("bcs.jpeg")
    bcs_plot = ggplot(met_data, aes(x=bcs, y=g2_measure)) + 
    geom_violin() + 
    geom_point(size = 3) + 
    xlab("BCS") + 
    annotate(geom = 'text', label = text_p, x = 0.75, y = 17.8, hjust = 0, vjust = 1, cex = 8) + 
    ggtitle("B") + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 23, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05), 
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 23, color = "black", face = "bold")) + 
    stat_summary(fun=mean, geom="crossbar", width = 0.2, color = "#0072B2")
    print(bcs_plot)
dev.off()

# weight
results = lm(met_data$wt2~met_data$g2_measure)
s = summary(results)
my.p = s$coefficients[2,4]
my.p = round(my.p, digits = 3)
text_p=paste(toString("p="), toString(my.p), sep="")
jpeg("weight.jpeg")
    weight_plot = ggplot(met_data, aes(x= (wt2), y=g2_measure)) + 
    ggtitle("C") + 
    xlab("Weight (kg)") + 
    annotate(geom = 'text', label = text_p, x = 454, y = 17.8, hjust = 0, vjust = 1, cex = 8) + 
    geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
    geom_point(size = 3) + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 23, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=22),
    axis.title.x = element_text(size = 23, color = "black", face = "bold")) + # also made the font one size smaller
    xlim(c(450,910)) # had to make the upper x limit greater than 2000 bc the 2000 was getting cut off of the graph. By making the upper x limit 2050 it extended the graph enough to not have it cut off :) 
    print(weight_plot)
dev.off()

# days_pp_fa
results = lm(met_data$days_pp_fa~met_data$g2_measure)
s = summary(results)
my.p = s$coefficients[2,4]
my.p = round(my.p, 3)
text_p=paste(toString("p="), toString(my.p), sep="")
jpeg("days_pp_fa.jpeg")
    days_pp_fa_plot = ggplot(met_data, aes(x=days_pp_fa, y=g2_measure)) + 
    ggtitle("D") + 
    xlab("Days Postpartum") + 
    ylab("Follicle Size (mm)") + 
    annotate(geom = 'text', label = text_p, x = 51, y = 17.8, hjust = 0, vjust = 1, cex = 8) + 
    geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
    geom_point(size = 3) + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 23, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05),
    axis.title = element_text(size = 23, color = "black", face = "bold")) 
    print(days_pp_fa_plot)
dev.off()

# pg_to_g2

results = lm(met_data$pg_to_g2~met_data$g2_measure)
s = summary(results)
my.p = s$coefficients[2,4]
my.p = round(my.p, 3)
text_p=paste(toString("p="), toString(my.p), sep="")
jpeg("pg_to_g2.jpeg")
    pg_to_g2_plot = ggplot(met_data, aes(x=pg_to_g2, y=g2_measure)) + 
    ggtitle("E") + 
    xlab("PGF to GnRH2 (hrs)") + 
    annotate(geom = 'text', label = text_p, x = 48.6, y = 17.8, hjust = 0, vjust = 1, cex = 8) + 
    geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
    geom_point(size = 3) + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 23, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05), 
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 23, color = "black", face = "bold")) 
    print(pg_to_g2_plot)
dev.off()



# g2_to_fa

results = lm(met_data$g2_to_fa~met_data$g2_measure)
s = summary(results)
my.p = s$coefficients[2,4]
my.p = round(my.p, 3)
text_p=paste(toString("p="), toString(my.p), sep="")
jpeg("g2_to_fa.jpeg")
    g2_to_fa_plot = ggplot(met_data, aes(x=g2_to_fa, y=g2_measure)) + 
    ggtitle("F") + 
    xlab("GnRH2 to Aspiration (hrs)") + 
    annotate(geom = 'text', label = text_p, x = 16.38, y = 17.8, hjust = 0, vjust = 1, cex = 8) + 
    geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
    geom_point(size = 3) + 
    theme_bw() +
    theme(axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 23, color = "black"),
    axis.text = element_text(size = 23, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    plot.title = element_text(hjust = -0.05), 
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 23, color = "black", face = "bold")) 
    print(g2_to_fa_plot)
dev.off()

# combining all of the figures into a gridded figure

jpeg(filename = 'cow_traits.jpeg', width = 1800, height = 1200)
grid.arrange(age_plot, bcs_plot, weight_plot, days_pp_fa_plot, pg_to_g2_plot, g2_to_fa_plot, nrow = 2)
dev.off()