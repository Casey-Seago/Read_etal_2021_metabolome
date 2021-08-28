###################################################################################
############################ Heatmap Figure Generation ############################
###################################################################################

# load packages
library(ggplot2)
library(sqldf)
library(grid)
library(gridExtra)  # to arrange individual plots into a grid of figures
library(ComplexHeatmap)
library(RColorBrewer)
# import data 
metab = read.csv(file="metab_cleaned.csv")
fdr = read.csv(file = "all_pvalues.csv")

# creating a heatmap figure of all metabolites
# removing the first column that is unimportant and that R adds when saved as a csv
metab = metab[,2:40]

# transposing the dataframe
tmetab = t(metab)

# creating a version where the dataframe is put in acsending follicle size order
tometab = tmetab[,order(tmetab[1,])]

# turning the row that is g2_measure into the column names so that the follicle size is the label for the column
tometab = data.frame(tometab, stringsAsFactors = FALSE) # needs to have stringsAsFactors = FALSE
colnames(tometab) <- tometab[1,] # setting the column names to the values in the first column
tometab <- tometab[-1, ] # deleting the first row bc it is the g2_measure

# ordering the metabolites in order of fdr value and then by alphabetical order of metabolite (some have the same fdr value and it randomly assigns orders to them. Want to keep it consistent between figures)
fdr$metabolite = as.character(fdr$metabolite)
ordered_fdr = order(fdr$fdrtest, fdr$metabolite) # decreasing = true bc of how I want them to be ordered on the final heatmap

# creating a column separated string that I can use to reorder the rows
ordered_fdr_string = paste(ordered_fdr, collapse = ",")

# reordering the dataframe
tometab = tometab[c(15,17,18,9,13,4,33,20,24,1,26,35,5,22,11,23,32,6,7,36,16,29,2,28,21,12,30,3,19,27,8,10,38,37,14,34,25,31),]

# creating a character vector for improved row name appearance 
fn_or = function(x) gsub("\\.", "/", x) # . is a wildcard symbol, have to escape that with the \\
fn_hyph = function(x) gsub("\\.", "-", x)
fn_rem = function(x) gsub("X", "", x)
fn_acid = function(x) gsub("-acid", " acid", x)
metab_names = list() # creating an empty list for correctly formatted names to be added to
rownames(tometab) = as.factor(rownames(tometab)) # converting to factor bc needed classification for following code
for (i in 1:nrow(tometab)) {
    if (startsWith(rownames(tometab)[i], "xylose")) {
    } else if (startsWith(rownames(tometab)[i], "X")) {
        metab = fn_acid(fn_hyph(fn_rem(rownames(tometab)[i])))
    } else if (startsWith(rownames(tometab)[i], "valine" ) | 
        startsWith(rownames(tometab)[i], "leucine") | 
        startsWith(rownames(tometab)[i], "citrate") | 
        startsWith(rownames(tometab)[i], "homoserine") | 
        startsWith(rownames(tometab)[i], "alanine")) {
        metab = fn_acid(fn_or(rownames(tometab)[i]))
    } else {
        metab = fn_acid(fn_hyph(rownames(tometab)[i]))
    }
    metab_names = append(metab_names, metab)
}

# replacing the 0 values with NA prior to log2() transformation of data
tometab[ tometab == 0] <-NA
tometab = log2(tometab)

# converting to a matrix for the heatmap function
tometabM = data.matrix(tometab)

jpeg('heatmap_ordered_scaled.jpeg', height = 1800, width = 1800)
 Heatmap(tometabM, name = "log2(peak intensity value)", 
    col = brewer.pal(9,"YlOrRd"), 
    row_labels = metab_names, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    row_names_side = "left", 
    width = unit(40, "cm"), height = unit(50, "cm"), 
    heatmap_legend_param = list(title_position = "lefttop-rot"), 
    column_names_gp = grid::gpar(fontsize = 23),
    row_names_gp = grid::gpar(fontsize = 23), 
    na_col = "white")
dev.off()
