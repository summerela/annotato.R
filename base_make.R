#!/usr/bin/env Rscript
############################################################################################
###############################  base_make.R  ###############################################
##  PURPOSE: This program takes as input a fasta and the fastq files it was generated from.
##           It runs the programs needed to generate a bam alignment file, protein
##           translations and blast reports needed for downstream analysis and annotation.
##  Output:
##       bowtie.bam = sorted bowtie alignment file
##       bowtie.csorted.bam = coordinant sorted bam filex
##       fasta_blastx = blastx report of nucleotide fasta against protein db, tab separated 
##       RSEM.genes.results = RSEM gene level abundance estimates
##       RSEM.isoforms.results = RSEM isoform level abundance estimates
##
##       The final output file will be placed in output_dir/base_files/
############################################################################################
options(error=traceback)
##########################
## Read in config file  ##
#########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

##check that output directory exists
dir.create(file.path(output_dir))

##create base_files directory
dir.create(file.path(paste(output_dir, "base_files", sep="/")))
base_output <- paste(output_dir, "base_files", sep="/")

##change working directory to output directory
setwd(base_output)

#################
##RSEM command ##
#################
##RSEM paired-end reads
if (reads == "PAIRED"| reads == "paired"){
    ##if fastqs are zipped, unzip
    if (length(grep("*.gz$", left)) > 0){
        system(paste("parallel gunzip :::", left, right, sep=" "))
        left <- lapply(strsplit(left, ".gz$"), function(x) x[1])
        right <- lapply(strsplit(right, ".gz$"), function(x) x[1])
    }
    ##check for strandedness
    if (direction == "1"){
        rsem_cmd <- paste(rsem_path, "--transcripts", final.fa, "--seqType fq --left", left,
                          "--right", right, "--est_method RSEM --aln_method bowtie --SS_lib_type RF --thread_count",
                          CPU, "--output_dir", base_output, "--coordsort_bam --prep_reference --trinity_mode", sep=" ")
    }
    if (direction == "0"){
        rsem_cmd <- paste(rsem_path, "--transcripts", final.fa, "--seqType fq --left", left,
                          "--right", right, "--est_method RSEM --aln_method bowtie --SS_lib_type FR --thread_count",
                          CPU, "--output_dir", base_output, "--coordsort_bam --prep_reference --trinity_mode", sep=" ")
    }
    if (direction == "2"){
        rsem_cmd <- paste(rsem_path, "--transcripts", final.fa, "--seqType fq --left", left, "--right",
                          right, "--est_method RSEM --aln_method bowtie --thread_count", CPU,
                          "--output_dir", base_output, "--coordsort_bam --prep_reference --trinity_mode", sep=" ")
    }
}
##RSEM single end reads
if (reads == "SINGLE" | reads == "single"){
    ##check for compressed fq files
    if (length(grep("*.gz$", single)) > 0){
        system(paste("gunzip", single, sep=" "))
        single <- lapply(strsplit(single, ".gz$"), function(x) x[1])
    }
    if (direction == "1"){
        rsem_cmd <- paste(rsem_path, "--transcripts", final.fa, "--seqType fq --single", single,
                          "--est_method RSEM --aln_method bowtie --SS_lib_type R --thread_count", CPU,
                          "--output_dir", base_output, "--coordsort_bam --prep_reference --trinity_mode", sep=" ")
    }
    if (direction == "0"){
        rsem_cmd <- paste(rsem_path, "--transcripts", final.fa, "--seqType fq --single", single,
                          "--est_method RSEM --aln_method bowtie --SS_lib_type F --thread_count", CPU,
                          "--output_dir", base_output, "--coordsort_bam --prep_reference --trinity_mode", sep=" ")
    }
    if (direction == "2"){
        rsem_cmd <- paste(rsem_path, "--transcripts", final.fa, "--seqType fq --single",
                          single, "--est_method RSEM --aln_method bowtie --thread_count", CPU,
                          "--output_dir", base_output, "--coordsort_bam --prep_reference --trinity_mode", sep=" ")
    }
}
##############################################
## Filter out discordant and unmapped reads ##
##############################################
sam_filter <- function(bam_file){
  ##remove unmapped reads
  system(paste("samtools view -b -F 4", bam_file, ">", paste(base_output, "final.fa_mapped.bam", sep="/")), wait=TRUE)
  ##keep only properly paired reads
  system(paste("samtools view -b -f 0x02", paste(base_output, "final.fa_mapped.bam", sep="/"), ">", paste(base_output, "final.fa_paired.bam", sep="/")), wait=TRUE)
  ##remove multi-mapped reads by skipping alignments with mapq < 1
  system(paste("samtools view -bq 15", paste(base_output, "final.fa_paired.bam", sep="/"), ">", paste(base_output, "final.fa_final.bam", sep="/")), wait=TRUE)
  ##index file for viewing in IGV 
  system(paste("samtools index", paste(base_output, "final.fa_final.bam", sep="/")), wait=TRUE) 
  ##create counts file
  system(paste("samtools idxstats", paste(base_output, "final.fa_final.bam", sep="/"), ">", paste(base_output, "final.fa_counts.txt", sep="/")), wait=TRUE)
}

####################
## Blastx command ##
####################
blastx_cmd <-  paste("blastx -db", blastdb, "-query", final.fa, "-outfmt 6 -num_threads", CPU, "-evalue", eval, "-max_target_seqs 1 -out", paste(base_output, paste(basename(final.fa), "blastx", sep="_"), sep="/"))

########################################
## Blastx against chicken for TransPS ##
########################################
blastx_chick<-  paste("blastx -db /home/sre/BlastDB/chickenDB -query", final.fa, "-outfmt 6 -num_threads", CPU, "-evalue", eval, "-max_target_seqs 1 -out", paste(base_output, paste(basename(final.fa), "blastx_chick", sep="_"), sep="/"))

#########################
## create qsub command ##
#########################
qsub_cmd <- paste("qsub -N Base_Files -j y -pe by_node_32", CPU, "-p -100 -l mem_free=4G,h_vmem=4.2G -M", email, "-m abe -b y", sep=" ")

#####################
## submit command ##
#####################
if (qsub == "true" | qsub == "TRUE"){
    cat("\n Submitting RSEM...\n")
    system(paste(qsub_cmd, rsem_cmd, sep=" "))
    cat("\n Filtering out unmapped and discordant reads...\n")
    sam_filter(paste(qsub_cmd, paste(base_output, "bowtie.csorted.bam", sep="/")))
    cat("\n Submitting blastx...\n")
    system(paste(qsub_cmd, blastx_cmd, sep=" "))
    cat("\n Submitting blastx_chicken...\n")
    system(paste(qsub_cmd, blastx_chick, sep=" "))
    cat("Finished submitting jobs...\n")
}
if (qsub == "false" | qsub == "FALSE"){
    cat("\n Running RSEM...\n")
    system(paste(rsem_cmd), wait=TRUE)
    cat("\n Filtering out unmapped and discordant reads...\n")
    sam_filter(paste(base_output, "bowtie.csorted.bam", sep="/"))
    cat("\n Running blastx...\n")
    system(paste(blastx_cmd))
    cat("\n Running blastx_chicken...\n")
    system(paste(blastx_chick))
    cat("Finished running jobs. Please check stdout for job status messages...\n")
}




