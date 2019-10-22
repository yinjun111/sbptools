#

#R script to filter file for PCA

library("argparser",quietly =T)

version="0.1"


description=paste0("count2cpm\nversion ",version,"\n","Usage:\nDescription: Convert Count to CPM\n")


#####
#Add arguments
#####

parser <- arg_parser(description=description)

parser <- add_argument(parser, arg="--count", type="character", help = "Count file")
parser <- add_argument(parser, arg="--cpm", type="character", help = "CPM file")


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


convert_count2cpm <- function(mat,na.rm=0) {
	#remove NA
	mat[is.na(mat)]<-na.rm
	
	mat.sum<-apply(mat,2,sum)

	#counts per million
	mat.cpm<-t(apply(mat,1,function(x) {x/(mat.sum/1000000)}))
	
	return(mat.cpm)
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


#filter based on count
data.cpm<-convert_count2cpm(data.count)

#
write_table_proper(file=args$cpm,data.cpm,"Feature")

