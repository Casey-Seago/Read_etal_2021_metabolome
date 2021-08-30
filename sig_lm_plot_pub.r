###################################################################################################
############################ Significant Metabolite LM Plot Generation ############################
###################################################################################################

# load packages
library(ggplot2)
library(sqldf)
library(grid)
library(gridExtra)  # to arrange individual plots into a grid of figures
library(tidyr)
library(stringr)


# import data as metab
metab = read.csv(file = "metab_cleaned.csv")
fdr_sig = read.csv(file = "fdr_sig_pvalues.csv")

# generating lm plots of significant metabolites with the pvalue on each chart

# removing the . from the metabolite names in the fdr_sig dataframe as well as the metab dataframe bc sql does not like . being in the column names
fdr_sig$metabolite = gsub("\\.","_",fdr_sig$metabolite)
colnames(metab) = gsub("\\.","_",colnames(metab))
# creating a new dataframe with only columns for the significant metabolites and g2_measure
sig_metab = sqldf(paste("SELECT", "g2_measure,", paste(fdr_sig$metabolite, collapse = ","), "FROM metab", sep = " "))
sig_metab = droplevels(sig_metab)

# formatting the metabolite names
# have names that need to be metabolite_metabolite -> metabolite/metabolite, Xmetabolite -> metabolite, and metabolite_metabolite to metabolite-metabolite
fn_or = function(x) sub("_", "/", x)
fn_hyph = function(x) gsub("_", "-", x)
fn_rem = function(x) sub("X", "", x)
fn_acid = function(x) sub("-acid", " acid", x)

# reordering the dataframe here so that the panel labels get placed on the correct plot
ordered_fdr = order(fdr_sig$fdrtest,fdr_sig$metabolite) # ordering by significance and, in the event fdr is the same, by alphabetical order
metab_df_order = ordered_fdr + 1 # did plus one bc the metab dataframe has the first column as g2_measure
metab_df_order_string = paste(metab_df_order, collapse = ",")
sig_metab = sig_metab[,c(1,9,10,11,6,8,3,18,12,15,2,16,19,4,13,7,14,17,5)] # ordering the df based on significance and alphabetical order

#reordering the FDR reference df
ordered_fdr_string = paste(ordered_fdr, collapse = ",")
fdr_sig = fdr_sig[c(8,9,10,5,7,2,17,11,14,1,15,18,3,12,6,13,16,4),]

# generating the individual plots
plot_list = list() # creating an empty list for plot names to be added to
for (i in 2:ncol(sig_metab)) {
    my_fdr = round(fdr_sig$fdrtest[i-1], digits=3)
    text_fdr=paste(toString("FDR="), toString(my_fdr), sep="")
    if (startsWith(colnames(sig_metab[i]), "X")){
        metab = fn_rem(colnames(sig_metab)[i])
        metab = fn_hyph(metab)
        metab = fn_acid(metab)
    } else if (startsWith(colnames(sig_metab[i]), "valine")) {
        metab = fn_or(colnames(sig_metab)[i])
        metab = fn_acid(metab)
    } else if (startsWith(colnames(sig_metab[i]), "leucine")) {
        metab = fn_or(colnames(sig_metab)[i])
        metab = fn_acid(metab)
    } else {
        metab = fn_hyph(colnames(sig_metab)[i])
        metab = fn_acid(metab)
    }
    plot_name = paste(colnames(sig_metab[i]), "plot", sep ="_")
    plot_list = append(plot_list, plot_name) # generating a list of plot names to use later during grid.arrange()
    panel = LETTERS[i-1] # using LETTERS function in r to get the panel labelling for the figure. Have to subtract one bc you want to start with A, but your column index is 2.
    assign(plot_name, ggplot(sig_metab, aes_string(x= colnames(sig_metab[1]), y=colnames(sig_metab[i]))) + ####!!! Have to have aes_string() instead of aes() bc of the way you are referencing your column names. I put the constant, x = g2_measure column in the same format as the changing y = colnames(sig_metab[i]) for ease/consistency
        ggtitle(panel) +
        xlab("") +
        ylab("") +
        annotate(geom = 'text', label = metab, x = Inf, y = -Inf, hjust = 1, vjust = -1.8, cex = 6.5, fontface = 2) +
        annotate(geom = 'text', label = text_fdr, x = Inf, y = -Inf, hjust = 1, vjust = -0.3, cex = 6.5) +
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
        axis.text.x = element_text(size=22)))
}

# creating a comma separated string for use in grid.arrange()
plot_string = paste(plot_list, collapse = ",")

# arranging the plots into a grid of figures and outputting as a jpeg
jpeg(filename = 'metabolite_lm_plots.jpeg', width = 1800, height = 2600)
    grid.arrange(alpha_ketoglutarate_plot,glutamate_plot,O_acetyl_L_serine_plot,N_acetyl_beta_alanine_plot,aspartate_plot,creatinine_plot,D_gluconate_plot,methionine_plot,phenylalanine_plot,pyruvate_plot,uric_acid_plot,uridine_plot,X2_oxoisovalerate_plot,X3_methylphenylacetic_acid_plot,leucine_isoleucine_plot,allantoin_plot,X2_dehydro_D_gluconate_plot,valine_betaine_plot, 
        nrow = 6, 
        top = "", 
        bottom = textGrob("Follicle Diameter at GnRH2 (mm)", 
            gp=gpar(fontsize=24, font =2)), 
        left = textGrob("Peak Intensity", 
            rot = 90, 
            gp=gpar(fontsize=24, font = 2)))
dev.off()
