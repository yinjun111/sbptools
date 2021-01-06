version <- 0.21

#v0.21, add pdf for the plots

# ---------------------
# Required libraries
# ---------------------
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("ggfortify"))
suppressPackageStartupMessages(library("Cairo"))
suppressPackageStartupMessages(library("reshape2"))

# ---------------------
# Options parsing
# ---------------------
option_list <- list(
  make_option(c("--input"), default="gene.results.merged.expr.filtered.auto.txt",
              help="The expression matrix to be analyzed. [default=%default]"),
  make_option(c("--config"), default="Config.txt",
              help="Configuration file with group information. [default=%default]"), 
  make_option(c("--group"), default="Group",
              help="Group name. [default=%default]"),   
  make_option(c("--geneanno"), default="geneanno.txt",
              help="Gene annotation file. [default=%default]"),   
  make_option(c("--pn"), default="RNA_",
              help="Project name. [default=%default]"),  
  make_option(c("--out"), default="EXPR-QC",
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

parser <- OptionParser(usage = "%prog [options] expr", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

# ---------------------
# Output directory
# ---------------------
# Create output directory
dir.create(file.path(opt$out), showWarnings = T)

# ---------------------
# Log file
# ---------------------
logfile <- paste0(opt$out, "/", opt$pn, "_expr-qc_run-log")
sink(logfile)

# Add file and options log 
cat(paste0("Start time: ", Sys.time(), "\n\n"))
cat(paste0("Current version: [", version, "]\n\n"))
cat("Following parameters used:\n")
cat(paste0("    Project name: [", opt$pn,"]\n"))
cat(paste0("    Output directory: [", opt$out, "]\n"))
cat(paste0("    Expression matrix converted to Log2: [", opt$log2, "]\n"))
cat(paste0("    Psuedocount used for Log2: [", opt$pseudocount, "]\n"))
cat(paste0("    Group info: [", opt$group, "]\n"))
cat(paste0("    Expression matrix: [", opt$input, "]\n"))
cat(paste0("    Configuration file: [", opt$config, "]\n"))
cat(paste0("    Gene annotation: [", opt$geneanno, "]\n"))
cat(paste0("    Make boxplot for genes: [", opt$boxplot, "]\n\n"))
cat(paste0("expr-qc version ", version, " running ...\n\n"))

# ---------------------
# Input expr & Config
# ---------------------
expr <- read.delim(opt$input, sep = "\t", header = T, row.names = 1, check.names = F)
expr[is.na(expr)] <- 0
config <- read.delim(opt$config, sep = "\t", header = T)
geneanno <- read.delim(opt$geneanno, sep = "\t", header = T, stringsAsFactors = F)

# ---------------------------------
# Define groups, colors, Log2 norm
# ---------------------------------
# Define group information for samples
Group <- factor(config[,opt$group][match(colnames(expr), config$Sample)])
pca.groups <- data.frame(Group)

cat("Following samples found:\n")
cat(do.call(paste, c(as.list(config$Sample), sep = ", ")))
cat("\n\n")
cat("Following groups found:\n")
cat(do.call(paste, c(as.list(levels(Group)), sep = ", ")))
cat("\n\n")

# Define colors for plots (up to 10 sample groups)
plot_cols <- c("#EA3323", "#0E2BF5", "#41ab5d", "#c994c7","#962219", 
               "#A3FCFE", "#A1F96F", "#D373DA", "#ED712D", "#F4B86B")
cols <- plot_cols[1:nlevels(Group)]

# Convert to Log2 if TRUE
if (opt$log2 == T) {
  expr.pca <- expr + opt$pseudocount
  expr.pca <- log2(expr.pca)
}

# ---------------------
# Correlation plot
# ---------------------
cat("Creating correlation plot ...\n")
M <- cor(expr)
cor_plot <- paste0(opt$out,"/" , opt$pn, "_correlation-plot.png")
CairoPNG(filename = cor_plot,res = 300,width=2500, height=2200, pointsize = 140/nrow(config))
corrplot(M, method = "color", addCoef.col="grey", type = "upper")
dev.off()

# ---------------------
# Sample boxplot
# ---------------------
cat("Creating boxplots for all samples ...\n")
# Melt expr data frame and add group information
expr.melt <- melt(expr.pca)
expr.melt$Group <- config[, opt$group][(match(expr.melt$variable, config[, "Sample"]))]

p <- ggplot(expr.melt, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=Group, alpha=0.75), outlier.shape=NA)  + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_fill_manual(values=cols) + xlab("Sample") + ylab("Expression level (Log2)") + labs(title="") +
  scale_alpha(guide = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
box_plot <- paste0(opt$out,"/" , opt$pn, "_sample-boxplots_", ".png")
CairoPNG(filename = box_plot,res = 300,width=1000 + 500*dim(config)[1]/5, height=2200)
print(p)
dev.off()

# ---------------------
# PCA
# ---------------------
# Compute prcomp
expr.pca <- t(expr.pca)
expr.pca <- expr.pca[, apply(expr.pca, 2, var, na.rm = T) != 0]
expr.pca <- prcomp(expr.pca, center = T, scale. = T)

# Plot of variances associated with each PC
cat("Creating plot of variances associated with each PC ...\n")
var_plot <- paste0(opt$out,"/" , opt$pn, "_PCA-variances.png")
CairoPNG(filename = var_plot,res = 300,width=2500, height=2200)
plot(expr.pca, type = "l")
dev.off()

# PC1 vs PC2
cat("Creating PCA plot PC1 vs. PC2 with sample names ...\n")
pca_plot1 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&2-names.png")
CairoPNG(filename = pca_plot1,res = 300,width=2500, height=2200)
autoplot(expr.pca, x=1, y=2, data=pca.groups, colour="Group", label = TRUE, label.size = 4, shape=F) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols)
dev.off()


cat("Creating PCA plot PC1 vs. PC2 ...\n")
pca_plot2 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&2.png")
CairoPNG(filename = pca_plot2,res = 300,width=2500, height=2200)
autoplot(expr.pca, x=1, y=2, data=pca.groups, colour="Group", size = 4, alpha = 0.75) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 
dev.off()

pca_plot2.pdf <- paste0(opt$out,"/" , opt$pn, "_PCA-1&2.pdf")
pdf(file =pca_plot2.pdf,width=7.5, height=6.6) #default 7x7 inch
autoplot(expr.pca, x=1, y=2, data=pca.groups, colour="Group", size = 4, alpha = 0.75) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 
dev.off()


# PC1 vs PC3
cat("Creating PCA plot PC1 vs. PC3 with sample names ...\n")
pca_plot3 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&3-names.png")
CairoPNG(filename = pca_plot3,res = 300,width=2500, height=2200)
autoplot(expr.pca, x=1, y=3, data=pca.groups, colour="Group", label = TRUE, label.size = 4, shape=F) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols)
dev.off()

cat("Creating PCA plot PC1 vs. PC3 ...\n")
pca_plot3 <- paste0(opt$out,"/" , opt$pn, "_PCA-1&3.png")
CairoPNG(filename = pca_plot3,res = 300,width=2500, height=2200)
autoplot(expr.pca, x=1, y=3, data=pca.groups, colour="Group", size = 4, alpha = 0.75) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 
dev.off()

# ---------------------
# Boxplots of genes
# ---------------------
cat("Creating boxplots for individual genes ...\n")
genes <- unlist(strsplit(opt$boxplot, ","))
if (length(genes) > 0){
  for (gene in genes){
    gene.ensembl <- geneanno$Gene[match(gene, geneanno$gene_name)]
    gene.expr <- unlist(expr[gene.ensembl, ])
    gene.df <- data.frame(expr=gene.expr, group=config[, opt$group])
    
    p <- ggplot(gene.df, aes(x=group, y=expr)) + 
      geom_boxplot(aes(fill=group, alpha=0.75), outlier.shape=NA)  + theme_classic() + 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="grey28", alpha=1)+
      scale_fill_manual(values=cols) + xlab("Group") + ylab("Expression level (Log2)") + labs(title=gene) +
      scale_alpha(guide = 'none') + theme(plot.title = element_text(size=22))
    box_plot <- paste0(opt$out,"/" , opt$pn, "_Boxplot_",gene, ".png")
    CairoPNG(filename = box_plot,res = 300,width=2500, height=2200)
    print(p)
    dev.off()
  }
}