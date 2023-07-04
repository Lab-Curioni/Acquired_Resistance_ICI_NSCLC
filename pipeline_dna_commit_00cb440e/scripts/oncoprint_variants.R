args <- commandArgs(trailingOnly = TRUE)
# oncoprint of variant in samples, 
# seiler_tcga_2017
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")

library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)

variantInfo_file= args[1]
outFile_filteredTable = args[2]
outFile_oncoprint = args[3]
filterThreshold = as.integer(args[4])
isDrugRelated = args[5] #yes/no
inputFile_mutLoad_all = args[6]
excludeTTN = args[7] #yes/no

#variantInfo_file = "/Users/fzickman/Desktop/projects/seiler/tcga/results/variantInfo/variantInfo_allSamples_drug_related_impact26.tsv"
variantTable = read.table(variantInfo_file,header = TRUE,stringsAsFactors=FALSE, sep = "\t")
dfVariants = data.frame(variantTable,check.names = FALSE)
dim(dfVariants) #2420

dfVariants_filtered = subset(dfVariants, Num_mutated_samples >= filterThreshold)
dim(dfVariants_filtered)

#excludeTTN = "yes"
if (excludeTTN == "yes"){
	dfVariants_filtered = subset(dfVariants_filtered, dfVariants_filtered$Gene != "TTN")
	dim(dfVariants_filtered)
}

# write out filtered table
#outTable_filtered =  paste(outFile_filteredTable,sep="")
write.table(dfVariants_filtered,file=outFile_filteredTable,sep="\t", row.names = FALSE, quote = FALSE)

# prepare for oncoprint
dfVariants_filtered[is.na(dfVariants_filtered)] = ""
#names(dfVariants_filtered)
dfVar_mat = cbind(dfVariants_filtered$Gene,dfVariants_filtered[,3:length(names(dfVariants_filtered))])

dfVar_mat[dfVar_mat=="0"]<-""

rownames(dfVar_mat) = dfVar_mat[, 1]  # to remove first column, first set the names as rownames
dfVar_mat = dfVar_mat[, -1]  # now remove column

# convert the 5 and 3 prime to "Five" and "Three" prime
dfVar_mat_converted <- t(apply(dfVar_mat, 1, function(this_row){ 
  new_row <- vapply(this_row, function(y){
    these_entries <- unlist(strsplit(y, ";"))
    new_entries <- c()
    for(this_entry in these_entries) {
      if(grepl("prime", this_entry)){
        this_entry <- sub("5_", "five_", sub("3_", "three_", this_entry))
      } 
      new_entries <- c(new_entries, this_entry)
    }
    paste0(new_entries, collapse=";")
  }, character(1))
  return(new_row)
}))

# colors and maps
these_cols_raw <- c(brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"), 
                    brewer.pal(12, "Paired"), brewer.pal(8, "Set2"))

#barplot(rep(1, length(these_cols_raw)), col = these_cols_raw, names=seq(1, length(these_cols_raw)), cex.names = 0.5)

these_cols <- these_cols_raw[c(26, 8, 32, 6, 24, 18, 22, 20, 34, 7, 17, 33, 23, 21, 25, 19, 27)]
#barplot(rep(1, length(these_cols)), col = these_cols, names=c("Start_lost", "Stop_lost", "Stop_gained", "Start_gained", "Frameshift_indel",
             #                                                 "Inframe_indel", "Missense", "Splice_site", "Protein_interacion_loci", "5_prime_UTR", "3_prime_UTR", "Synonymous",
             #                                                 "Intronic", "Nc_region", "Sequence_feature"), las = 2, cex.names = 0.8)

cols_mutation_type = c(start_lost = these_cols[1],
                       stop_retained_variant = these_cols[2],
                       stop_lost = these_cols[6], 
                       stop_gained = these_cols[8],
                       start_gained = these_cols[5], 
                       frameshift_variant = these_cols[3],
                       #Inframe_indel = these_cols[6],
                       missense_variant = these_cols[7],
                       splice_donor_variant = these_cols[4],
                       splice_acceptor_variant = these_cols[9],
                       splice_region_variant = these_cols[14],
                       #Protein_interacion_loci = these_cols[9],
                       five_prime_UTR_variant = these_cols[10],
                       three_prime_UTR_variant = these_cols[11], 
                       synonymous_variant = these_cols[12],
                       intron_variant = these_cols[13],
                       #Nc_region = these_cols[14],
                       disruptive_inframe_deletion = these_cols[15],
                       conservative_inframe_deletion = these_cols[16],
                       disruptive_inframe_insertion = these_cols[17],
                       MULTIPLE = "black")
scale_radius_MULTIPLE <- 0.3

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  ####################################
  start_lost = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["start_lost"], col = NA))
  },
  stop_lost = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["stop_lost"], col = NA))
  },
  stop_gained = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["stop_gained"], col = NA))
  },
  start_gained = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["start_gained"], col = NA))
  },
  frameshift_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["frameshift_variant"], col = NA))
  },
  missense_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["missense_variant"], col = NA))
  },
  splice_donor_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["splice_donor_variant"], col = NA))
  },
  splice_acceptor_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["splice_acceptor_variant"], col = NA))
  },
  splice_region_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["splice_region_variant"], col = NA))
  },
  five_prime_UTR_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["five_prime_UTR_premature_start_codon_gain_variant"], col = NA))
  },
  three_prime_UTR_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["three_prime_UTR_variant"], col = NA))
  },
  synonymous_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["synonymous_variant"], col = NA))
  },
  intron_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["intron_variant"], col = NA))
  },
  disruptive_inframe_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["disruptive_inframe_deletion"], col = NA))
  },
  conservative_inframe_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["conservative_inframe_deletion"], col = NA))
  },
  stop_retained_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["stop_retained_variant"], col = NA))
  },
  disruptive_inframe_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cols_mutation_type["disruptive_inframe_insertion"], col = NA))
  },
  ##
  # MULTIPLE
  MULTIPLE = function(x, y, w, h) {
    w = convertWidth(w, "cm")*0.9
    h = convertHeight(h, "cm")*0.9
    l = min(unit.c(w, h)) # l*0.5 is radius
    grid.circle(x, y, l*scale_radius_MULTIPLE, gp = gpar(fill = cols_mutation_type["MULTIPLE"], col = NA))
  }
)


titleString = "Frequently mutated genes"

if (isDrugRelated == "yes") {
  titleString = "Frequently mutated drug-related genes"
}


#outfile_oncoprint_png = paste(outFile_oncoprint,sep="")
#png(filename=outFile_oncoprint, width = 4000, height = 4500, res = 300)

#ht = oncoPrint(dfVar_mat_converted, get_type = function(x) strsplit(x, ";")[[1]],
#               alter_fun = alter_fun, col = cols_mutation_type, 
#               column_title = titleString,
#               remove_empty_columns = TRUE,
#               #show_column_names = TRUE,
#               barplot_ignore = "MULTIPLE",
#               #row_order = NULL, column_order = NULL,
#               gap = unit(3, "mm"),
#               row_names_gp = gpar(fontsize = 7),
#               row_barplot_width = unit(2, "cm"),
#               show_pct = TRUE,
#               heatmap_legend_param = list(title = "Alterations", at = c(names(cols_mutation_type)), 
#                                           labels = names(cols_mutation_type), nrow = 3)) #,
#nrow = 1, title_position = "leftcenter"),
#heatmap_legend_param = list(title = "Drug", at = c("present"), 
#                           labels = c("Present"), nrow = 1, title_position = "leftcenter"))

#draw(ht, heatmap_legend_side = "bottom")

#dev.off()

#inputFile_mutLoad_all = "/Users/fzickman/Desktop/projects/seiler/tcga/results/variantInfo/samples_variantInfo_all.txt"
varInfo_table_all = read.table(inputFile_mutLoad_all,header=TRUE,sep="\t",check.names = FALSE)
dfVarInfo_all = data.frame(varInfo_table_all,check.names = FALSE)

#dfSig = data.frame(patient = dfVarInfo_all$patient,signature = dfVarInfo_all$signature)
dfSig = data.frame(signature = dfVarInfo_all$signature)
head(dfSig)
rownames(dfSig) = colnames(dfVar_mat_converted)
head(dfSig)

#barplot(rep(1, length(these_cols_raw)), col = these_cols_raw, names=seq(1, length(these_cols_raw)), cex.names = 0.5)

sig_type_colors <-these_cols_raw[c(5,6,7,8,11,17,9,13,14,22,24,26,34)]  # 13 signatures
col_assign_sig_types <- setNames(sig_type_colors, unique(dfSig$signature))

#barplot(rep(1, length(col_assign_sample_types)), col = sample_type_colors, names=unique(dfSig_tcga$signature), las = 2, cex.names = 0.8)

png(filename=outFile_oncoprint, width = 4000, height = 4500, res = 300)

htSig = oncoPrint(dfVar_mat_converted, get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = cols_mutation_type, 
               column_title = titleString,
               remove_empty_columns = TRUE,
               #show_column_names = TRUE,
               barplot_ignore = "MULTIPLE",
               #row_order = NULL, column_order = NULL,
               gap = unit(3, "mm"),
               row_names_gp = gpar(fontsize = 7),
               row_barplot_width = unit(2, "cm"),
               show_pct = TRUE,
               heatmap_legend_param = list(title = "Alterations", at = c(names(cols_mutation_type)), 
                                           labels = names(cols_mutation_type), nrow = 3),
               bottom_annotation = HeatmapAnnotation(df = dfSig, 
                                                     col = list(signature = col_assign_sig_types)))

draw(htSig, heatmap_legend_side = "bottom")

dev.off()
