# sbptools

sbptools is a suite of computational pipelines to analyze NGS data.

It is free to use and distribute this software package, but without any warranty.

sbptools
version: $version
Usage: sbptools [tool] [parameters]


Parameters:

    ########
    #Supported in both Firefly and Falco
    ########
    rnaseq-process    RNA-seq QC, Align, and RSEM for FASTQ files
    rnaseq-merge      Merge rnaseq-process results for downstream analyses
    rnaseq-de         Perform DE analysis using DESeq2
    rnaseq-summary    Summarize RNA-Seq DE results

    rnaseq-var        RNA-seq variant calling pipeline
    rnaseq-motif      RNA-seq TFBS motif finding pipeline
    rnaseq-motif-summary  RNA-seq TFBS motif finding results summary	
    motif-finder      Transcription factor binding motif prediction
	
    chipseq-process   ChIP-seq/ATAC-Seq QC, Align, and Peak Calling
    chipseq-merge     Summarize ChIP-Seq results
    chipseq-de        Perform DE analysis for ChIP-Seq
    chipseq-summary   Summarize ChIP-seq DE results

    dnaseq-process    DNA-Seq (Exome/Genome-Seq) processing based on GATK4

    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file (by Andrew Hodges)

    gsea-gen          Generate files ready for GSEA analysis

    ########
    #Supported only in Firefly
    ########
    parallel-job      Batch job submission in Firefly

    ########
    #Supported only in Falco
    ########	
	
    bs-fastq          Download and merge FASTQ files from Basespace	
    ensembl2ucsc      Convert Ensembl gtf/fasta/bed into UCSC format
    geo-download      Download raw FASTQ files from GEO
