#!/opt/patric-common/runtime/bin/Rscript

#parameter format: 
#RunDESeq2.R <counts_file.txt> <metadata_file.txt> <output_prefix> <htseq|stringtiet> <contrast 1> <contrast 2> ... <contrast n>
#contrasts should be a csv pair and list all the contrasts to make in this dataset
args = commandArgs(trailingOnly=TRUE)

# TODO: adjust for removing report img path name

numContrasts = length(args) - 5

if (numContrasts < 1) {
    stop("Not enough parameters: run_deseq2.R <counts_file> <metadata_file> <output_prefix> <report_img_path> <bacteria|host> <contrast 1> ... <contrast n>")
}
counts.file = args[1]
metadata.file = args[2]
output_prefix = args[3]
report_prefix = args[4]
genome_type = args[5]

#Check file extensions
if ("bacteria" == genome_type) {
    count_sep = "\t"
} else if ("host" == genome_type) {
    count_sep = "," 
} else { 
    print("Error in run_deseq2.R: can't determine counts file delimeter")
    print(counts.file)
    stop() 
}

#Differential expression library, load quietly so it doesn't completely fill the log output
suppressMessages(library(DESeq2,quietly=TRUE))
#library(ggplot2,quietly=TRUE)
#library(EnhancedVolcano,quietly=TRUE)
#library(gridExtra,quietly=TRUE)
#library(svglite)

#Load counts table and metadata table 
count.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
#rownames(count.mtx) = gsub("gene-","",rownames(count.mtx))
#rownames(count.mtx) = gsub("rna-","",rownames(count.mtx))
metadata <- read.table(metadata.file,sep='\t',header=T,stringsAsFactors=TRUE)

print(head(count.mtx))
print(head(metadata))
print(metadata$Condition)
print(metadata$Sample)

#calculate png width and height
#png_width = 600
#if (numContrasts > 1) {
#    png_width = png_width * 2
#}
#png_height = 460
#if (numContrasts > 1) {
#    png_height = ceiling((numContrasts/2)) * png_height
#}
#svg_width = 14
#svg_height = ceiling((numContrasts/2)) * 5

#iterate over contrasts
#index 5 in args is where the contrasts currently start
plot_list = vector("list",numContrasts)
contrast.index <- 0
for (i in 6:length(args)) 
{
    contrast.index <- contrast.index + 1
    #Subset data on current contrast
    curr_contrast = unlist(strsplit(args[i],":"))
    #curr_contrast = gsub("-","_",curr_contrast)
    curr.metadata = subset(metadata,(subset=Condition==curr_contrast[1])|(subset=Condition==curr_contrast[2]))
    curr.count.mtx = count.mtx[,curr.metadata$Sample] 
    #TODO(wrong?): Remove all zero rows and add pseudocount to genes: Will not work for some samples otherwise
    curr.count.mtx = curr.count.mtx[rowSums(curr.count.mtx) > 0,]
    #curr.count.mtx = curr.count.mtx + 1
    #Run standard DESeq2 pipeline
    print(paste("running DESeq: ",curr_contrast[1]," against ",curr_contrast[2],sep=""))
    dds <- DESeqDataSetFromMatrix(countData = curr.count.mtx, colData = curr.metadata, design = ~Condition)
    dds <- DESeq(dds)
    res <- results(dds,contrast=c("Condition",curr_contrast[1],curr_contrast[2]))

    res = cbind(res,data.frame(Gene_Name=rownames(res)))
    res = res[,c("Gene_Name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    #write to output file
    results_file = paste(output_prefix,curr_contrast[1],"_vs_",curr_contrast[2],".deseq2.tsv",sep="")     
    write.table(res,file=results_file,sep="\t",quote=FALSE,row.names=FALSE)

    #Create volcano plot
    # moved to different file
    #min_x_axis = min(res$log2FoldChange) - 1
    #max_x_axis = max(res$log2FoldChange) + 1
    #ev_image_name = paste(output_prefix,"_",curr_contrast[1],"_vs_",curr_contrast[2],"_mqc.png",sep="")
    #contrast_name = paste(curr_contrast[1]," over ",curr_contrast[2],sep="")
    #png(ev_image_name,width=png_width,height=png_height)
    #ev_img <- EnhancedVolcano(res,lab=rownames(res),x='log2FoldChange',y='padj',xlim=c(min_x_axis,max_x_axis),subtitle="",title=contrast_name,legendPosition="top",titleLabSize=14)
    #plot_list[[contrast.index]] <- ev_img
    #print(ev_img)
    #dev.off()
} 

###Output PNG
#grid_png = paste("Volcano_Plots_mqc.png",sep="")
#png(grid_png,width=png_width,height=png_height)
#do.call("grid.arrange",c(plot_list,ncol=2))
#dev.off()

###Output SVG
#grid_svg = paste("Volcano_Plots_mqc.svg",sep="")
#grid_svg = paste(report_prefix,"Volcano_Plots.svg",sep="")
#svg(grid_svg,width=svg_width,height=svg_height)
#svglite(grid_svg,width=svg_width,height=svg_height)
#do.call("grid.arrange",c(plot_list,ncol=2))
#dev.off()
