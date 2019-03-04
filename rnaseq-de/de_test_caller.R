
library(DESeq2)



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

deseq2_test <- function(mat,anno,design,fc_cutoff=1,q_cutoff=0.05,pmethod="Wald",qmethod="BH",core=4){
  
  #anno and design may need to be checked


	dds <- DESeqDataSetFromMatrix(countData = round(mat),
	                              colData = anno,
	                              design= as.formula(design))
	dds <- DESeq(dds,test=pmethod)
	resultsNames(dds) # lists the coefficients
	res <- results(dds, pAdjustMethod =qmethod)
	
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


setwd("D:\\Work\\SBP\\sbptools-dev")

data<-read.table("gene.results.merged.count.txt",header=T,row.names=1,sep="\t")
anno<-read.table("Sample_Anno.txt",header=T,row.names=1,sep="\t")


data.sel<-filter_data(data)

data.sel.result<-deseq2_test(data.sel,anno=anno,design="~Tissue")
