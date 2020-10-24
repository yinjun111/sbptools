#!/usr/bin/env Rscript

# ---------------------
# Required libraries
# ---------------------
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("ggfortify"))
suppressPackageStartupMessages(library("Cairo"))

# ---------------------
# Options parsing
# ---------------------
option_list <- list(
  make_option(c("--input"), default="gene.results.merged.tpm.filtered.auto.txt",
              help="The TPM matrix to be analyzed. [default=%default]"),
  make_option(c("--config"), default="Config.txt",
              help="Configuration file with group information. [default=%default]"), 
  make_option(c("--geneanno"), default="geneanno.txt",
              help="Gene annotation file. [default=%default]"),   
  make_option(c("--pn"), default="RNA_",
              help="Project name. [default=%default]"),  
  make_option(c("--out"), default="TPM-QC",
              help="Output folder name. [default=%default]"),
  make_option(c("--log2"), action="store_true", default=TRUE, 
              help="Apply log2. NAs are treated as 0s and a pseudocount is added if specified. [default=%default]"),
  make_option(c("--pseudocount"), type="double", default=1e-02,
              help="Specify a pseudocount for log2 conversion. [default=%default]"),
  make_option(c("--boxplot"), default="",
              help="Make boxplots for comma separated list of genes. [default=%default]"),  
  make_option(c("-v", "--verbose"), default=FALSE, action="store_true",
              help="Verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] TPM", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

# Create output directory
dir.create(file.path(opt$out), showWarnings = T)

# ---------------------
# Input TPM & Config
# ---------------------
tpm <- read.delim(opt$input, sep = "\t", header = T, row.names = 1)
tpm[is.na(tpm)] <- 0

config <- read.delim(opt$config, sep = "\t", header = T)
geneanno <- read.delim(opt$geneanno, sep = "\t", header = T, stringsAsFactors = F)

# ---------------------
# Correlation plot
# ---------------------
M <- cor(tpm)
cor_plot <- paste0(opt$out,"/" , opt$pn, "_correlation-plot.png")
CairoPNG(filename = cor_plot,res = 300,width=2500, height=2200, pointsize = 7)
corrplot(M, method = "color", addCoef.col="grey", type = "upper")
dev.off()

# ---------------------
# 2D PCA
# ---------------------
if (opt$log2 == T) {
  tpm.pca <- tpm + opt$pseudocount
  tpm.pca <- log2(tpm.pca)
}
tpm.pca <- t(tpm.pca)
tpm.pca <- tpm.pca[, apply(tpm.pca, 2, var, na.rm = T) != 0]
tpm.pca <- prcomp(tpm.pca, center = T, scale. = T)

# Plot of variances associated with each PC
var_plot <- paste0(opt$out,"/" , opt$pn, "_PCA-variances.png")
CairoPNG(filename = var_plot,res = 300,width=2500, height=2200)
plot(tpm.pca, type = "l")
dev.off()

# Define group information for samples
Group <- factor(config$Group[match(rownames(tpm.pca$x), config$Sample)])
pca.groups <- data.frame(Group)

# Define colors for plots (up to 10 sample groups)
plot_cols <- c("#EA3323", "#0E2BF5", "#41ab5d", "#c994c7","#962219", 
               "#A3FCFE", "#A1F96F", "#D373DA", "#ED712D", "#F4B86B")
cols <- plot_cols[1:nlevels(Group)]

# PC1 vs PC2
pca_plot1 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&2-names.png")
CairoPNG(filename = pca_plot1,res = 300,width=2500, height=2200)
autoplot(tpm.pca, x=1, y=2, data=pca.groups, colour="Group", label = TRUE, label.size = 4, shape=F) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols)
dev.off()

pca_plot2 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&2.png")
CairoPNG(filename = pca_plot2,res = 300,width=2500, height=2200)
autoplot(tpm.pca, x=1, y=2, data=pca.groups, colour="Group", size = 4, alpha = 0.75) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 
dev.off()

# PC1 vs PC3
pca_plot3 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&3-names.png")
CairoPNG(filename = pca_plot3,res = 300,width=2500, height=2200)
autoplot(tpm.pca, x=1, y=3, data=pca.groups, colour="Group", label = TRUE, label.size = 4, shape=F) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols)
dev.off()

pca_plot3 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&3.png")
CairoPNG(filename = pca_plot3,res = 300,width=2500, height=2200)
autoplot(tpm.pca, x=1, y=3, data=pca.groups, colour="Group", size = 4, alpha = 0.75) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 
dev.off()

# Make boxplots
genes <- unlist(strsplit(opt$boxplot, ","))
if (length(genes) > 0){
  for (gene in genes){
    gene.ensembl <- geneanno$Gene[match(gene, geneanno$gene_name)]
    #gene.tpm <- log2(unlist(tpm[gene.ensembl, ]) + opt$pseudocount)
    gene.tpm <- unlist(tpm[gene.ensembl, ])
    gene.df <- data.frame(tpm=gene.tpm, group=config$Group)
    
    p <- ggplot(gene.df, aes(x=group, y=tpm)) + 
      geom_boxplot(aes(fill=group, alpha=0.75), outlier.shape=NA)  + theme_classic() + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="black", alpha=1)+
      scale_fill_manual(values=cols) + xlab("Group") + ylab("TPM") + labs(title=gene) +
      scale_alpha(guide = 'none') + theme(plot.title = element_text(size=22))
    
    box_plot <- paste0(opt$out,"/" , opt$pn, "_Boxplot_",gene, ".png")
    CairoPNG(filename = box_plot,res = 300,width=2500, height=2200)
    print(p)
    dev.off()
  }
}
