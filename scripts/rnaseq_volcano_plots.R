#!/opt/patric-common/runtime/bin/Rscript

#parameter format:
#rnaseq_volcano_plots.R <output_prefix> <deseq2_file1> <contrast_name1> <deseq2_file2> <contrast_name2>... 
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
    stop("Incorrect number of arguments: rnaseq_volcano_plots.R <output_prefix> <deseq2_file1> <contrast_name1> <deseq2_file2> <contrast_name2>...")
}

library(gridExtra,quietly=T)
library(ggplot2,quietly=T)
library(EnhancedVolcano,quietly=T)
library(svglite)

output_prefix = args[1]
output_name = paste(output_prefix,'volcano_plot.svg',sep='')
arg_pairs = args[2:length(args)]
plot_index <- 0
deseq_files <- arg_pairs[c(TRUE,FALSE)]
plot_list = vector("list",length(deseq_files))
print(deseq_files)
contrast_names <- arg_pairs[c(FALSE,TRUE)]
print(contrast_names)
curr_index <- 0
# If file exists assume it's the correct file
for (diffexp_file in deseq_files) {
    curr_index <- curr_index + 1
    if (file.exists(diffexp_file)) {
        plot_index <- plot_index + 1
        curr.data <- read.table(diffexp_file,sep='\t',header=T)
        contrast_title <- contrast_names[curr_index]
        ev_img <- EnhancedVolcano(curr.data,lab=curr.data$Gene_Name,x='log2FoldChange',y='padj',subtitle="",title=contrast_title,legendPosition="top",titleLabSize=14) 
        plot_list[[plot_index]] <- ev_img
    }
}

svg_width = 10
numContrasts = length(deseq_files)
svg_height = numContrasts * 5
svglite(output_name,width=svg_width,height=svg_height)
do.call("grid.arrange",c(plot_list,ncol=1))
dev.off()
