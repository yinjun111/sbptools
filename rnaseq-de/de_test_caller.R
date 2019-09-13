

library("argparser",quietly =T)

version="0.3"

#0.2b, change auto filter to *5. Add indfilter and cookscutoff option
#0.23, add write_table_proper
#0.3, add MA plot,fixed xlim for volcano plot

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
parser <- add_argument(parser, arg="--independentfiltering",type="logical", help = "DESeq2 independentFiltering",default=T)
parser <- add_argument(parser, arg="--cookscutoff",type="logical", help = "DESeq2 cooksCutoff",default=T)


args = parse_args(parser)

print(args)

#other parameters

core=1

#args_file = "tempArgObjectFile.rds"
#saveRDS(args, args_file); print(args); quit(); #comment this after creating args_file
#args = readRDS(args_file)  

#independentfiltering=T
#cookscutoff=T


#####
#Filter genes
#####


filter_data <- function(mat,type="sum",cutoff=10,na.rm=0) {
	#remove NA
	mat[is.na(mat)]<-na.rm


	if(cutoff=="auto") {
	  num.cutoff<-ncol(mat)*5 # count cutoff eq two times number of samples #changed to *5 3/26
	} else {
		if(grepl("^auto",cutoff,perl=T)) {
			fold<-regmatches(cutoff,regexpr("\\d+$",cutoff,perl=T))
			num.cutoff<-ncol(mat)*as.numeric(fold)
		} else {
			num.cutoff<-as.numeric(cutoff)
		}
	}

		
	if(type == "sum") {
		mat.sel<-mat[apply(mat,1,sum)>=num.cutoff,]
	} else if (type=="none") {
		mat.sel<-mat
	}
	
	return(mat.sel)
}


######
#DESeq2
######

deseq2_test <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,pmethod="Wald",qmethod="BH",core=core,treat,ref,independentfiltering=T,cookscutoff=T){
  library(DESeq2,quietly =T)
  #anno and design may need to be checked


	dds <- DESeqDataSetFromMatrix(countData = round(mat),
	                              colData = anno,
	                              design= as.formula(design))
	dds <- DESeq(dds,test=pmethod)
	resultsNames(dds) # lists the coefficients
	
	comp<-tail(all.vars(as.formula(design)),n=1) #the last variable of the formula is used as comparison
	
	if(cookscutoff) {
		res <- results(dds, pAdjustMethod =qmethod,contrast=c(comp,treat,ref),independentFiltering=independentfiltering)
	} else {
		res <- results(dds, pAdjustMethod =qmethod,contrast=c(comp,treat,ref),independentFiltering=independentfiltering,cooksCutoff=cookscutoff)	
	}
	
	fc<-res[,2]
	p<-res[,5]
	q<-res[,6]
	stat<-apply(res[,c(1,3,4)],1,function(x) {paste(x,collapse = ",")})
	
	#significance by fc & qval
	sig<-rep(0,length(q))
	
	sig[!is.na(fc) & fc>=fc_cutoff & !is.na(q) & q<q_cutoff]=1
	sig[!is.na(fc) & fc<=-fc_cutoff & !is.na(q) & q<q_cutoff]=-1
	sig[is.na(fc) | is.na(q) ]=NA

	mat.result<-cbind(fc,stat,p,q,sig)
	
	colnames(mat.result)<-c(values(res)[[2]][2],"DESeq2 Stat:Mean,SE,Wald stat",values(res)[[2]][5],values(res)[[2]][6],paste("Significance: Log2FC ",round(fc_cutoff,3)," ",qmethod, "P ",q_cutoff,sep=""))
	
	
	
	result<-list()
	result$result=mat.result
	result$dds=dds
	result$res=res

	return(result)
}

volcano_plot_ggplot<-function(fc,q,sig,xlim=c(-5,5),ylim=c(0,20),xlab="Log2FC",ylab="-log10 P",main="Volcano Plot",fc_cutoff=args$fccutoff,q_cutoff=args$q_cutof) {
  
  library(ggplot2)
  
  q[is.na(q)]<-1 #remove NA
  
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
  
  
  #redefine cols and shs
  cols.sel<-cols[levels(as.factor(sig.new))]
  shs.sel<-shs[levels(as.factor(shapes))]
  
  #print(cols.sel)
  #print(shs.sel)
  
  #transform data
  fc[fc>xlim[2]]=xlim[2]
  fc[fc<xlim[1]]=xlim[1]
  
  q[-log10(q)>ylim[2]]=10^-ylim[2]
  
  #defined cols and shapes
  
  data<-data.frame( lfc=fc,q=-log10(q),sig=sig.new,shape=shapes)
  
  
  #plot
  vol <- ggplot(data, aes(x = lfc, y =q, fill = sig ,shape=factor(shape)))
  
  vol + ggtitle(label = main) +
    geom_point(size = 2, alpha = 1, na.rm = T, colour = "black") +
    scale_fill_manual(name="Color",values = cols) +
    scale_shape_manual(name="Shape",values = shs) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) + # change the legend
    guides(fill=guide_legend(override.aes = list(size = 3,  colour=cols.sel),title="Significance"),shape="none")+
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    ylim(ylim) +
    scale_x_continuous(breaks =seq(xlim[1],xlim[2],1),lim=xlim) +
    geom_hline(yintercept = -log10(as.numeric(q_cutoff)), colour="#990000", linetype="dashed") + #p cutoff
    geom_vline(xintercept = as.numeric(fc_cutoff), colour="#990000", linetype="dashed") + geom_vline(xintercept = -as.numeric(fc_cutoff), colour="#990000", linetype="dashed") + # fc cutoff line
    annotate("text",x=xlim[2]-0.5, y=ylim[2]-0.5, label=length(which(sig==1)),colour = "red",size = 5) + 
    annotate("text",x=xlim[1]+0.5, y=ylim[2]-0.5, label=length(which(sig==-1)),colour = "green",size = 5) + 
    annotate("text",x=xlim[2]-1, y=-log10(as.numeric(q_cutoff))+0.5, label=substring(main,28),colour = "red",size = 3) 
  
}


#to be implemented
ma_plot_ggplot<-function(m,a,sig,xlim=c(1,20),ylim=c(-5,5),xlab="A:Log2 Mean of Normalized Counts",ylab="M:Log2 Fold Change",main="MA Plot",m_cutoff=args$fc_cutof) {
  
  #m is lfc
  #a is mean
  
  library(ggplot2)
  
  #a[is.na(a)]<-1 #remove NA
  
  m<-as.numeric(unlist(m))
  a<-as.numeric(unlist(a))
  
  #define color  
  cols <- c("Up" = "red", "Down" = "green","N.S."="grey")
  shs <- c("21" = 21, "24" = 24)
  #define col and shape
  
  shapes=rep(21,length(m))
  shapes[abs(m)>ylim[2]]=24
  shapes[a>xlim[2]]=24
  shapes[a<xlim[1]]=24
  
  sig.new<-rep("N.S.",length(sig))
  sig.new[sig==1]="Up"
  sig.new[sig==-1]="Down"
  
  
  #redefine cols and shs
  cols.sel<-cols[levels(as.factor(sig.new))]
  shs.sel<-shs[levels(as.factor(shapes))]
  
  #print(cols.sel)
  #print(shs.sel)
  
  #transform data
  a[a>xlim[2]]=xlim[2]
  a[a<xlim[1]]=xlim[1]
  
  m[m>ylim[2]]=ylim[2]
  m[m<ylim[1]]=ylim[1]
  
  #defined cols and shapes
  
  data<-data.frame( m=m,a=a,sig=sig.new,shape=shapes)
  
  
  #plot
  vol <- ggplot(data, aes(x = a, y =m, fill = sig ,shape=factor(shape)))
  
  vol + ggtitle(label = main) +
    geom_point(size = 2, alpha = 1, na.rm = T, colour = "black") +
    scale_fill_manual(name="Color",values = cols) +
    scale_shape_manual(name="Shape",values = shs) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) + # change the legend
    guides(fill=guide_legend(override.aes = list(size = 3,  colour=cols.sel),title="Significance"),shape="none")+
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    ylim(ylim) +
    scale_x_continuous(breaks =seq(xlim[1],xlim[2],1),lim=xlim) +
    #geom_hline(yintercept = -log10(as.numeric(a_cutoff)), colour="#990000", linetype="dashed") + #p cutoff
    geom_hline(yintercept = as.numeric(m_cutoff), colour="#990000", linetype="dashed") + geom_hline(yintercept = -as.numeric(m_cutoff), colour="#990000", linetype="dashed") + # a cutoff line
    annotate("text",x=xlim[2]-0.5, y=ylim[2]-0.5, label=length(which(sig==1)),colour = "red",size = 5) + 
    annotate("text",x=xlim[2]-0.5, y=ylim[1]+0.5, label=length(which(sig==-1)),colour = "green",size = 5) 
  
}

write_table_proper<-function(file,data,name="Gene") {
	data.df<-data.frame(name=rownames(data),data)
	names(data.df)<-c(name,colnames(data))

	write.table(data.df,file=file, row.names=FALSE,sep="\t",quote=F)

	#write.table(data.frame(name=rownames(data),data),file=file, row.names=FALSE,sep="\t",quote=F)
}

#generate plots

#X 1. MA, 
#2. volcano
#3. hist for fc?
#4. ...


#####
#Statistical test
#####


data<-read.table(args$"in",header=T,row.names=1,sep="\t",check.names=F,flush=T)
anno<-read.table(args$anno,header=T,row.names=1,sep="\t",check.names=F,flush=T)

#can be customized later
data.sel<-filter_data(data,cutoff=args$filter)

#save image
rdatafile=sub(".txt$",".rdata",args$out,perl=T)


data.sel.result<-deseq2_test(mat=data.sel,anno=anno,design=args$formula,fc_cutoff=args$fccutoff,q_cutoff=args$qcutoff,pmethod=args$pmethod,qmethod=args$qmethod,treat=args$treat,ref=args$ref,independentfiltering=args$independentfiltering,cookscutoff=args$cookscutoff)

#write.table(data.sel.result$result,file=args$out,sep="\t",quote=F,col.names = NA)

write_table_proper(data.sel.result$result,file=args$out,"Feature")

save.image(file=rdatafile)


#####
#Plots
#####

#volcano plot
vp_outfile_pdf=sub("\\.\\w+$","_volcanoplot.pdf",args$out,perl=T)

pdf(vp_outfile_pdf,width=7, height=6)
volcano_plot_ggplot(fc=data.sel.result$result[,1],q=data.sel.result$result[,4],sig=data.sel.result$result[,5],xlab=colnames(data.sel.result$result)[1],ylab=paste("-Log10 ",colnames(data.sel.result$result)[4],sep=""),main=paste("Volcano Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),q_cutoff=args$qcutoff,fc_cutoff = args$fccutoff)
dev.off()

#turned off for now, no X11 in linux

#vp_outfile_jpg=sub("\\.\\w+$","_volcanoplot.jpg",args$out,perl=T)

#jpeg(vp_outfile_jpg,width=7, height=6, units="in", res=300) 
#volcano_plot_ggplot(fc=data.sel.result$result[,1],q=data.sel.result$result[,4],sig=data.sel.result$result[,5],xlab=colnames(data.sel.result$result)[1],ylab=paste("-Log10 ",colnames(data.sel.result$result)[4],sep=""),main=paste("Volcano Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),q_cutoff=args$qcutoff,fc_cutoff = args$fccutoff)
#dev.off()


#####
#Plots
#####

#volcano plot
mp_outfile_pdf=sub("\\.\\w+$","_maplot.pdf",args$out,perl=T)

pdf(mp_outfile_pdf,width=7, height=6)
ma_plot_ggplot(m=data.sel.result$result[,1],a=log2(data.sel.result$res[,1]),sig=data.sel.result$result[,5],ylab=paste("M:",colnames(data.sel.result$result)[1],sep=""),main=paste("MA Plot ","Significance: Log2FC ",round(args$fccutoff,2)," ",args$qmethod, "P ",args$qcutoff,sep=""),m_cutoff = args$fccutoff)
dev.off()


#mp_outfile_jpg=sub("\\.\\w+$","_maplot.jpg",args$out,perl=T)

#jpeg(width=7, height=6, units="in", res=300)


