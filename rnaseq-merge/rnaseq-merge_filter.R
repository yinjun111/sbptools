#

#R script to filter file for PCA

library("argparser",quietly =T)

version="0.11"

#v0.11, add write_table_proper

description=paste0("rnaseq-merge_filter\nversion ",version,"\n","Usage:\nDescription: Filter merged results for PCA\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--count", type="character", help = "Count file")
parser <- add_argument(parser, arg="--fpkm", type="character", help = "FPKM file")
parser <- add_argument(parser, arg="--tpm", type="character", help = "TPM file")

#parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Annotation file")
parser <- add_argument(parser, arg="--filter",type="character", help = "Count filter for all the exps",default="auto")


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

write_table_proper<-function(file,data,name="Gene") {
	data.df<-data.frame(name=rownames(data),data)
	names(data.df)<-c(name,colnames(data))

	write.table(data.df,file=file, row.names=FALSE,sep="\t",quote=F)

	#write.table(data.frame(name=rownames(data),data),file=file, row.names=FALSE,sep="\t",quote=F)
}



#####
#Read files
#####

data.count<-read.table(args$count,header=T,row.names=1,sep="\t",check.names=F,flush=T)
data.fpkm<-read.table(args$fpkm,header=T,row.names=1,sep="\t",check.names=F,flush=T)
data.tpm<-read.table(args$tpm,header=T,row.names=1,sep="\t",check.names=F,flush=T)


#filter based on count
data.count.sel<-filter_data(data.count,cutoff=args$filter)
data.fpkm.sel<-data.fpkm[rownames(data.count.sel),]
data.tpm.sel<-data.tpm[rownames(data.count.sel),]

#
write_table_proper(file=sub(".txt",paste(".filtered.",args$filter,".txt",sep=""),args$count),data=data.count.sel,"Gene")
write_table_proper(file=sub(".txt",paste(".filtered.",args$filter,".txt",sep=""),args$fpkm),data=data.fpkm.sel,"Gene")
write_table_proper(file=sub(".txt",paste(".filtered.",args$filter,".txt",sep=""),args$tpm),data=data.tpm.sel,"Gene")

