#############################################################################
############################## Normality check ##############################
#############################################################################


# import data
data = read.csv(file="metab_cleaned.csv")

# removing the first column that R created when I saved it
data = data[,c(2:40)]

# checking the normality of the residuals by creating a qqplot. Outputting as a pdf so can just scroll through them
pdf('qq_metabolites.pdf')
for(i in 2:ncol(data)) { 
    txt = colnames(data[i])   
    model = lm(g2_measure ~ data[,i], data = data)
    plot.new()  
    ols_plot_resid_qq(model)
    text(x=.1, y=1, txt)
}
dev.off()

# doing a test for normality that will give us a shapiro wilk statistic and pvalue
for(i in 2:ncol(data)) {      
    sink(file = "resid_normality_check.txt", append = TRUE, type = c("output"), split = TRUE)
    model = lm(g2_measure ~ data[,i], data = data)
    print(colnames(data[i]))
    print(ols_test_normality(model))
    print('------------------------------------------------------------------')
    print('------------------------------------------------------------------')
    sink()
}

# the residuals are all approximately normal so no transformations are necessary