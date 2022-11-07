#!/opt/patric-common/runtime/bin/Rscript

#load libraries quietly
library(ggplot2,quietly=TRUE)
library(gridExtra,quietly=TRUE)
library(reshape2,quietly=TRUE)
library(svglite)

###Functions
#https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

###TODO: rename variables to not include "system"

#parameter format and check parameter inputs
#grid_violin_plots.R <system_map.csv> <counts_file.txt|csv> <output.file>
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Incorrect parameters: grid_violin_plots.R <gene_map.txt> <counts_file.txt|csv> <metadata.txt> <output_file_prefix>")
}

mapping.file = args[1]
counts.file = args[2]
metadata.file = args[3]
output.file = args[4]

#check files and extensions
#if (grepl("htseq",feature.count)) {
#    count_sep = "\t"
#} else if (grepl("stringtie",feature.count)) {
#    count_sep = ","
#} else {
#    print("Error in grid_violin_plots.R: can't determine counts file delimeter")
#    print(counts.file)
#    stop()
#}
#TPM matrix (replacing counts matrix) only uses tab delimeter at the moment
count_sep = "\t"

#read tables
counts.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
print('counts 1')
print(head(counts.mtx))
metadata <- read.table(metadata.file,sep="\t",header=T,stringsAsFactors=FALSE)
system.map <- read.table(mapping.file,sep="\t",header=T,stringsAsFactors=FALSE)

#Filter entries with no system label and get the intersection of patric_ids
system.map = system.map[!grepl("NONE",system.map[,2]),]
counts.mtx = counts.mtx[rownames(counts.mtx) %in% system.map[,1],]
print(head(system.map))
print('counts 2')
print(head(counts.mtx))

#Testing: min and max values
#log_min = log(min(counts.mtx)+1)
#log_max = log(max(counts.mtx))
min_val = min(counts.mtx)
max_val = max(counts.mtx)

#Get all unique systems and conditions
conditions = unique(metadata$Condition)
systems = unique(system.map[,2])  

#Calculate image width and height
#num_columns <- ceiling(sqrt(length(systems)))
num_columns <- 4 
num_samples <- ncol(counts.mtx)
num_rows <- ceiling(length(systems)/num_columns)
png_width = (num_columns + num_samples)*100
png_height = num_rows*200 
svg_width = num_columns + num_samples
svg_height = num_rows + 5

#create each plot and ad dit to a list of plots: do not render at this step: occurs when calling svglite() and do.call()
legend <- NULL 
plot_list = vector("list",length(systems)+1)
for (i in 1:length(systems)) {
    curr.system = systems[i] 
    curr.mtx = counts.mtx[rownames(counts.mtx) %in% system.map[which(system.map[,2] == curr.system),1],] 
    curr.mtx = data.frame(curr.mtx)
    curr.mtx$Genes <- rownames(curr.mtx)
    melt.df = melt(curr.mtx,id.vars=c("Genes"),measure.vars=colnames(curr.mtx)[-c(length(colnames(curr.mtx)))]) 
    colnames(melt.df) <- c("Gene","Sample","Counts")
    melt.df$LogCounts <- log(melt.df$Counts+1)
    melt.df$Condition <- rep(0,length.out=nrow(melt.df))
    for (c in conditions) {
        melt.df[melt.df$Sample %in% subset(metadata,subset=Condition==c)$Sample,]$Condition = c
    }
    x_axis_label = paste(toString(length(curr.mtx$Genes))," Genes",sep="")
    #vln_plot <- ggplot(melt.df,aes(x=Sample,y=LogCounts,fill=Condition))+geom_violin(trim=FALSE)+ylim(min_val,max_val)+ggtitle(curr.system)+ylab("TPM")+xlab(x_axis_label) 
    vln_plot <- ggplot(melt.df,aes(x=Sample,y=LogCounts,fill=Condition))+geom_violin(trim=FALSE)+ggtitle(curr.system)+ylab("TPM")+xlab(x_axis_label) 
    vln_plot = vln_plot + geom_boxplot(width=0.1,fill="white")
    if (is.null(legend)) {
        legend <- g_legend(vln_plot)
    }
    #vln_plot = vln_plot + theme(axis.text.x = element_text(size=8,angle=45,vjust=0.5), legend.position = "none")
    vln_plot = vln_plot + theme(axis.text.x = element_text(size=4,angle=315,vjust=0.5), legend.position = "none")
    plot_list[[i]] <- vln_plot
}
plot_list[[length(systems)+1]] <- legend

###Output PNG image
#vln_png = paste(output.file,"_Pathway_Distribution_mqc.png",sep="")
#png(vln_png,width=png_width,height=png_height)
#do.call("grid.arrange",c(plot_list,ncol=num_columns))
#dev.off()

#TODO: issue where it opens a second image and saves one as Rplot.pdf
###Output SVG image
#vln_svg = paste(output.file,"_Pathway_Distribution_mqc.svg",sep="")
vln_svg = paste(output.file,".svg",sep="")
#svg(vln_svg,width=svg_width,height=svg_height)
#svglite(vln_svg,width=svg_width,height=svg_height)
svglite(vln_svg,width=12,height=10)
do.call("grid.arrange",c(plot_list,ncol=num_columns))
dev.off()
