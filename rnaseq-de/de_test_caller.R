

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
	
	#need to modify here
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
	
	colnames(mat.result)<-c(values(res)[[2]][2],"DESeq2 Stat:Mean,SE,Wald stat",values(res)[[2]][5],values(res)[[2]][6],paste("Significance: Log2FC ",round(fc_cutoff,2)," ",qmethod, "P " ,q_cutoff,sep=""))
	
	return(mat.result)
}


#generate plots

#1. MA
#2. volcano
#3. hist for fc?
#4. ...


#####
#run
#####


data<-read.table(args$"in",header=T,row.names=1,sep="\t")
anno<-read.table(args$anno,header=T,row.names=1,sep="\t")


data.sel<-filter_data(data)

data.sel.result<-deseq2_test(mat=data.sel,anno=anno,design=args$formula,fc_cutoff=args$fccutoff,q_cutoff=args$q_cutoff,pmethod=args$pmethod,qmethod=args$qmethod,treat=args$treat,ref=args$ref)

write.table(data.sel.result,file=args$out,sep="\t",quote=F,col.names = NA)

