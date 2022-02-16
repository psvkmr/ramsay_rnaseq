### Analysis Pipeline: WS 


setwd("D:/Datasets/RNAseq/RNAseq_CD101_Tregs_WS")
fastq.files <- list.files('D:/Datasets/RNAseq/RNAseq_CD101_Tregs_WS/CD101_data-set/00_fastq', 
                          pattern="*.fastq.gz", 
                          full.names=TRUE)


####Quality control: FastQC####
#Quality control on the Fastq files has been performed in a previous pipeline
dir.create("fastqc")
cmd <- paste("fastqc --outdir fastqc", fastq.files)
for(file in seq(length(cmd))){
  print(cat("Analysis of:", cmd[file]))
  system(cmd[file])
  print("Analysis complete")
}


####Trimming: Trim Galore (wrapper for CutAdapt)####
for(idx in seq(1,length(fastq.files),2)){
  cmd <- paste("trim_galore --length 20 --paired --output_dir trimgalore", fastq.files[idx], fastq.files[idx + 1])
  system(cmd)
}


####Building reference transcriptome: Kallisto index####
# Creating index file from reference transcriptome (which will be used for pseudo alignment)
dir.create('Kallisto')
setwd("D:/Datasets/RNAseq/RNAseq_CD101_Tregs_WS/Kallisto")
system("wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")
# system("wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz") # non-coding RNA not used for building index
# system("cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.cDNA_ncRNA.fa.gz") # concatenate both files for full transcriptome. The "full" transcriptome was not use for the index

cmd <- "kallisto index -i GRCh38.99_cDNA_kallisto_index D:/Datasets/RNAseq/RNAseq_CD101_Tregs_WS/Kallisto/Reference_Genome/Homo_sapiens.GRCh38.cdna.all.fa.gz" # -o output directory not specified. will safe in current directory
system(cmd)



####Pseudo alignment: Mapping of trimmed reads to kallisto index####
# The output will be abundance.H5 files which can be imported via TxImport
fastq.files.trim <- list.files('/media/wlad/ILSE/RNAseq_data/trimgalore', pattern ='*.fq.gz',full.names=TRUE )
fastq.files.trim

fastq.names = sapply(strsplit(basename(fastq.files.trim), "\\."), '[[',1)
fastq.names = sapply(strsplit(fastq.names, "_"),'[[',1)
fastq.names
cmd <- sapply(seq(1,length(fastq.names),2), 
              function(idx)
                paste0('kallisto quant -i D:/Datasets/RNAseq/RNAseq_CD101_Tregs_WS/Kallisto/GRCh38.99_cnda_kallisto_index -o D:/Datasets/RNAseq/RNAseq_CD101_Tregs_WS/Kallisto/', 
                       fastq.names[idx] , " " , fastq.files.trim[idx] , " " , fastq.files.trim[idx+1]))
cmd
for(idx in cmd){
  system(idx)
}

