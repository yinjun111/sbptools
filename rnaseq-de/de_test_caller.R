

library("argparser",quietly =T)

version="0.1"

description=paste0("de_test\nversion ",version,"\n","Usage:\nDescription: Differential Expression calculation using DESeq2\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--in",short="-i", type="character", help = "Expr input file")
parser <- add_argument(parser, arg="--anno",short="-a", type="character", help = "Annotation file")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output file")
parser <- add_argument(parser, arg="--formula",short="-f", type="character", help = "DESeq formula")
parser <- add_argument(parser, arg="--treat",short="-t", type="character", help = "treatment name")
parser <- add_argument(parser, arg="--ref",short="-r", type="character", help = "reference name")
parser <- add_argument(parser, arg="--fccutoff", type="float", help = "Log2 FC cutoff",default=1)
parser <- add_argument(parser, arg="--qcutoff", type="float", help = "qcutoff",default=0.05)
parser <- add_argument(parser, arg="--pmethod",type="character", help = "Method used in DESeq2",default="Wald")
parser <- add_argument(parser, arg="--qmethod",type="character", help = "FDR Method",default="BH")
parser <- add_argument(parser, arg="--filter",type="float", help = "Count filter for all the exps",default="10")

args = parse_args(parser)

print(args)

#other parameters

core=1

#args_file = "tempArgObjectFile.rds"
#saveRDS(args, args_file); print(args); quit(); #comment this after creating args_file
#args = readRDS(args_file)  


#####
#Filter genes
#####


filter_data <- function(mat,type="sum",cutoff=10,na.rm=0) {
	#remove NA
	mat[is.na(mat)]<-na.rm
	
	if(type == "sum") {
		mat.sel<-mat[apply(mat,1,sum)>=cutoff,]
	} else if (type=="none") {
		mat.sel<-mat
	}
	
	return(mat.sel)
}


######
#DESeq2
######

deseq2_test <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,pmethod="Wald",qmethod="BH",core=core,treat,ref){
  library(DESeq2,quietly =T)
  #anno and design may need to be checked


	dds <- DESeqDataSetFromMatrix(countData = round(mat),
	                              colData = anno,
	                              design= as.formula(design))
	dds <- DESeq(dds,test=pmethod)
	resultsNames(dds) # lists the coefficients
	
	#need to modify here !!!
	comp<-substring(design,2)
	
	res <- results(dds, pAdjustMethod =qmethod,contrast=c(comp,treat,ref))
	
	fc<-res[,2]
	p<-res[,5]
	q<-res[,6]
	stat<-apply(res[,c(1,3,4)],1,function(x) {paste(x,collapse = ",")})
	
	#significance by fc & qval
	sig<-rep(0,length(q))
	
	sig[!is.na(fc) & fc>=fc_cutoff & !is.na(q) & q<q_cutoff]=1
	sig[!is.na(fc) & fc<=fc_cutoff & !is.na(q) & q<q_cutoff]=-1
	sig[is.na(fc) | is.na(q) ]=NA

	mat.result<-cbind(fc,stat,p,q,sig)
	
	colnames(mat.result)<-c(values(res)[[2]][2],"DESeq2 Stat:Mean,SE,Wald stat",values(res)[[2]][5],values(res)[[2]][6],paste("Significance: Log2FC ",round(fc_cutoff,2)," ",qmethod, "P ",q_cutoff,sep=""))
	
	return(mat.result)
}



volcano_plot_ggplot<-function(fc,q,sig,xlim=c(-5,5),ylim=c(0,20),xlab="Log2FC",ylab="-log10 P",main="Volcano Plot",fc_cutoff=args$fccutoff,q_cutoff=args$q_cutof) {

  fc<-as.numeric(unlist(fc))
  q<-as.numeric(unlist(q))
  
  #define color  
  cols <- c("Up" = "red", "Down" = "green","N.S."="grey")
  shs <- c("21" = 21, "24" = 24)
  #define col and shape
  
  shapes=rep(21,length(fc))
  shapes[abs(fc)>xlim[2]]=24
  shapes[-log10(q)>ylim[2]]=24
  
  sig.new<-rep("N.S.",length(sig))
  sig.new[sig==1]="Up"
  sig.new[sig==-1]="Down"
  
  #transform data
  fc[fc>xlim[2]]=xlim[2]
  fc[fc<xlim[1]]=xlim[1]
  
  p[-log10(p)>ylim[2]]=10^-ylim[2]
  
  #defined cols and shapes
  
  data<-data.frame( lfc=fc,q=-log10(q),sig=sig.new,shape=shapes)
  
  
  #plot
  vol <- ggplot(data, aes(x = lfc, y =q, fill = sig ,shape=factor(shape)))
  
  vol + ggtitle(label = main) +
    geom_point(size = 2, alpha = 1, na.rm = T, colour = "black") +
    scale_fill_manual(name="Color",values = cols) +
    scale_shape_manual(name="Shape",values = c(21,24)) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # change the legend
    guides(fill=guide_legend(override.aes = list(size = 3,  colour=c("green","grey","red"),title="Significance")),shape="none")+
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    ylim(ylim) +
    scale_x_continuous(breaks =seq(-5,5,1),lim=xlim) +
    geom_hline(yintercept = -log10(q_cutoff), colour="#990000", linetype="dashed") + #p cutoff
    geom_vline(xintercept = fc_cutoff, colour="#990000", linetype="dashed") + geom_vline(xintercept = -fc_cutoff, colour="#990000", linetype="dashed")  # fc cutoff line
    
}



#generate plots

#X 1. MA, 
#2. volcano
#3. hist for fc?
#4. ...


#####
#run
#####


data<-read.table(args$"in",header=T,row.names=1,sep="\t")
anno<-read.table(args$anno,header=T,row.names=1,sep="\t")

#need to be customized !!!
data.sel<-filter_data(data,cutoff=args$filter)

data.sel.result<-deseq2_test(mat=data.sel,anno=anno,design=args$formula,fc_cutoff=args$fccutoff,q_cutoff=args$q_cutoff,pmethod=args$pmethod,qmethod=args$qmethod,treat=args$treat,ref=args$ref)

write.table(data.sel.result,file=args$out,sep="\t",quote=F,col.names = NA)

#volcano plot
vp_outfile=sub("\\.\\w+$","_volcanoplot.pdf",args$out,perl=T)

pdf(vp_outfile)
volcano_plot_ggplot(fc=data.sel.result[,1],q=data.sel.result[,4],sig=data.sel.result[,5],xlab=colnames(data.sel.result)[1],ylab=paste("-Log10",colnames(data.sel.result)[4],sep=""),main="Volcano Plot",q_cutoff=args$q_cutoff,fc_cutoff = args$fccutoff)
dev.off()
