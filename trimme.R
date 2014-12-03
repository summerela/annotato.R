#!/usr/bin/env Rscript

##################################################################################
############################### trimme.R ######################################### 
##  PURPOSE: This script is used to trim fastq files before assembly. It takes as
##           input illumina reads and:
##                    - Removes illumina adapters for TruSeq3 (as used by HiSeq
#3                      and MiSeq machines) or  TruSeq2 (as used in GAII machines) 
##                    - Removes leading low quality or N bases (below quality 3)
##                    - Removes trailing low quality or N bases (below quality 3)
##                    - Scans the read with a 4-base wide sliding window,
##                      cutting when the average quality per base drops below 15
##                    - Drops reads below the 36 bases long
##                    - Trims the leading and trailing 5 base pairs
##
##  INPUT: single or paired end illumina fastq files specify options using config file
##  OUTPUT: trimmed reads, concatenated into one single-end file, or one left and one
##          right fastq file for paired-end reads
############################################################################################
##########################
## Read in config file  ##
##########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

##check if output directory exists, if not, create
dir.create(file.path(output_dir))

##create trimmed directory
dir.create(file.path(paste(output_dir, "trimmed", sep="/")))
trimmed_output <- paste(output_dir, "trimmed", sep="/")

##set working directory to output directory
setwd(output_dir)

##########################
## Setting up libraries ##
##########################
##check that packages are installed, and if not, install
packages <- c("parallel")
##load packages and supress their startup messages
suppressPackageStartupMessages(lapply(packages, require, character.only=T))

##########################
## Gather fastq files   ##
##########################
cat("Gathering fastq files... \n")
##data frame all files in directory path
fastq.df <- data.frame(path=list.files(fastq_dir, pattern="*.fastq.gz|fq.gz", full=TRUE))

##create list of adapters to pull
adapt_list <- gsub(",", "|", adapters)

##grap adapters from file names
fastq.df$adapter <- mclapply(strsplit(as.character(mclapply(strsplit(sub(".*/", "", fastq.df$path), "_"), function(x) x[4], mc.cores=5)), "[.]"), function(x) x[1], mc.cores=5)

##pull out fastqs that match adapters
fastq.df <- fastq.df[grep(adapt_list, fastq.df$adapter),]

##add left or right strand info
fastq.df$strand <- mclapply(strsplit(as.character(fastq.df$path), "_"), function(x) x[3], mc.cores=5)

##########################
## Merge fastq files   ##
##########################
if (reads == "paired"| reads == "PAIRED"){
    if (dim(fastq.df)[1] > 2){
        cat("\n Merging paired-end reads. \n")
        ##merge left reads
        system(paste("cat", paste(fastq.df[which(fastq.df$strand==1),]$path, collapse=" "), ">",
                     paste(trimmed_output, "merged_left.fq.gz", sep="/"), collapse=" "), wait=TRUE)
        cat("\n Uncompressing left reads...\n")
        ##gunzip left reads
        system(paste("gunzip", paste(trimmed_output, "merged_left.fq.gz", sep="/"), collapse=" "), wait=TRUE)
        ##merge right reads
        system(paste("cat", paste(fastq.df[which(fastq.df$strand==2),]$path, collapse=" "), ">", 
                     paste(trimmed_output, "merged_right.fq.gz", sep="/"), collapse=","), wait=TRUE)
        cat("\n Uncompressing right reads...\n")
        ##gunzip right reads
        system(paste("gunzip", paste(trimmed_output, "merged_right.fq.gz", sep="/"), collapse=" "), wait=TRUE)
    }
    if (dim(fastq.df)[1] <= 2){
        cat("\n Renaming paired-end reads. \n")
        ##merge left reads
        system(paste("cp", paste(fastq.df[which(fastq.df$strand==1),]$path, collapse=" "),
                     paste(trimmed_output, "merged_left.fq.gz", sep="/"), collapse=" "), wait=TRUE)
        cat("\n Uncompressing left reads...\n")
        ##gunzip reads
        system(paste("gunzip", paste(trimmed_output, "merged_left.fq.gz", sep="/"), collapse=" "), wait=TRUE)
        ##merge right reads
        system(paste("cp", paste(fastq.df[which(fastq.df$strand==2),]$path, collapse=" "),
                     paste(trimmed_output, "merged_right.fq.gz", sep="/"), collapse=","), wait=TRUE)
        cat("\n Uncompressing right reads...\n")
        ##gunzip right reads
        system(paste("gunzip", paste(trimmed_output, "merged_right.fq.gz", sep="/"), collapse=" "), wait=TRUE)
    }
}
if (reads == "single" | reads == "SINGLE"){
    if (dim(fastq.df)[1] >1){
        cat("\n Merging single-end reads. \n")
        ##merge left reads
        system(paste("cat", paste(fastq.df[which(fastq.df$strand==1),]$path, collapse=" "),
                     ">", paste(trimmed_output, "merged_single.fq.gz", sep="/"), collapse=" "), wait=TRUE)
        cat("\n Uncompressing reads...\n")
        ##gunzip reads
        system(paste("gunzip", paste(trimmed_output, "merged_single.fq.gz", sep="/"), collapse=" "), wait=TRUE)
    }
    if (dim(fastq.df)[1] == 1){
        cat("\n Renaming single-end read. \n")
        ##merge right reads
        system(paste("cp", paste(fastq.df[which(fastq.df$strand==2),]$path, collapse=" "),
                     paste(trimmed_output, "merged_single.fq.gz", sep="/"), collapse=","), wait=TRUE)
        cat("\n Uncompressing reads...\n")
        ##gunzip reads
        system(paste("gunzip", paste(trimmed_output, "merged_single.fq.gz", sep="/"), collapse=" "), wait=TRUE)
    }
}

##############################
## Create trimming command  ##
##############################

##set illumina adapter option
if (illumina == "2" & (reads == "paired" | reads == "PAIRED")){
    illumina.fa <- "/home/sre/tools/Trinity07172014/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic/adapters/TruSeq2-PE.fa"
}
if (illumina == "3" & (reads == "paired" | reads == "PAIRED")){
    illumina.fa <- "/home/sre/tools/Trinity07172014/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa"
}
if (illumina == "2" & (reads == "single" | reads == "SINGLE")){
    illumina.fa <- "/home/sre/tools/Trinity07172014/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic/adapters/TruSeq2-SE.fa"
}
if (illumina == "3" & (reads == "single" | reads == "SINGLE")){
    illumina.fa <- "/home/sre/tools/Trinity07172014/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic/adapters/TruSeq3-SE.fa"
}

##trim paired reads
if (reads == "paired" | reads == "PAIRED"){
    cat("\n Trimming paired-end reads...\n")
    trim <- paste("java -jar", trimmomatic, "PE -threads",
                  CPU, "-phred33 ./trimmed/merged_left.fq ./trimmed/merged_right.fq",
                  paste(trimmed_output, "trimmed_left.fq", sep="/"),
                  paste(trimmed_output, "left_unpaired.fq", sep="/"),
                  paste(trimmed_output, "trimmed_right.fq", sep="/"),
                  paste(trimmed_output, "right_unpaired.fq", sep="/"),
                  paste("ILLUMINACLIP:",
                        illumina.fa, ":2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:40:30 MINLEN:36",sep=""), sep=" ")
}

##trim single reads
if (reads == "single" | reads == "SINGLE"){
    cat("\n Trimming single-end reads...\n")
    trim <- paste("java -jar", trimmomatic, "SE -threads", CPU, "-phred33 ",
                  paste(trimmed_output, "merged_single.fq", sep="/"),
                  paste(trimmed_output, "trimmed_single.fq", sep="/"),
                  paste("ILLUMINACLIP:", illumina.fa,
                        ":2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36",sep=""), sep=" ")
}

#########################
## create qsub command ##
#########################
qsub_cmd <- paste("qsub -N Trim_Reads -j y -pe by_node_32", CPU, "-p -100 -l mem_free=4G,h_vmem=4.2G -M", email, "-m abe -b y", sep=" ")

#####################
## submit command ##
#####################
if (qsub=="true" | qsub == "TRUE"){
    cat("\n Running trimmomatic on the cluster...\n")
    system(paste(qsub_cmd, trim, sep=" "))
}
if (qsub == "false" | qsub == "FALSE"){
    cat("\n Running trimmomatic...\n")
    system(paste(trim))
}
           
cat("\n Trimming complete, files written to /trimmed folder.\n")



