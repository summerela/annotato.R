#!/usr/bin/env Rscript
#########################################################################################
###############################   stat_generato.R    ##################################### #########################################################################################
##  PURPOSE: This program takes as input a sorted bam file and fasta file,
##           and output reports on assembly quality.   
##
##  Output:
##       stats_fasta.fa = summary report of assembly, including:
##                           genes/transcript counts, GC%, N50 and length stats
##       fasta.fa_kmer_stats.txt = summary of kmer counts in assembly
##       kmer_histogram.jpg = graph of kmer count distribution
##       fasta.fa_length_dist.jpg = graph of assembly length distributiion
##       fasta.fa_blast_cov.w_pct_hit_length = report on blast hit coverage
##       blast_hit_cov.pdf = graph blast hit coverage of full and filtered assembly
##       bwt_align_stats.txt = bowtie alignment stats
#########################################################################################
options(error=traceback)
##########################
## Setting up libraries ##
##########################
##check that packages are installed, and if not, install
packages <- c(Biostrings")
##install any packages that are not currently installed
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.r-project.org")
}
##load packages
suppressPackageStartupMessages(lapply(packages, require, character.only=T))

#########################
## Read in config file  ##
##########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

##check that output directory exists
dir.create(file.path(output_dir))

##create base_files directory
dir.create(file.path(paste(output_dir, "qc", sep="/")))
qc_output <- paste(output_dir, "qc", sep="/")

##change working directory to output directory
setwd(qc_output)

##create of fasta files
fasta_files<- c(final.fa, filtered.fa)

##run full and filtered assembly through QC
run_QC <- function(fasta){
    ## run summary stats  ## 
    system(paste(fastatools, "stats", fasta, ">",paste(qc_output, paste("stats",
                                                                  basename(fasta), sep="_"), sep="/"), sep=" "))   
    ##Kmer Counts  ##
    ##create hash tables
    cat("\n \n Calculating kmer distribution... \n \n")
    system(paste("jellyfish count -m 25 -s 500000 -t", CPU, "-C", fasta, "-o",
                       paste("mer_counts", basename(fasta), sep="_"), sep=" "))  
    ##create list of hash tables 
    mer_counts <- list.files(qc_output, pattern=paste("mer_counts", basename(fasta), "*", sep="_"), full.names=TRUE)    
    ##merge hash tables if more than one
    if (length(mer_counts) > 1){
        system(paste("jellyfish merge -o counts.jf", paste(mer_counts, collapse=" "), sep=" "))
    }    
    else{
        system(paste("cp", mer_counts, paste(qc_output, paste(basename(fasta), "counts.jf", sep="_"), sep="/"), sep=" "))
    }
    ##count kmers
    system(paste("jellyfish stats",  paste(qc_output, paste(basename(fasta), "counts.jf", sep="_"), sep="/"), ">",
                 paste(qc_output, paste(basename(fasta), "kmer_stats.txt", sep="_"), sep="/"), sep=" "))
    ##create input for histogram
    system(paste("jellyfish histo -t", CPU, paste(qc_output, paste(basename(fasta), "counts.jf" ,sep="_"),
                                                  sep="/"), ">",
                 paste(qc_output, paste(basename(fasta), "kmer_hist.txt", sep="_"), sep="/"), sep=" "))
    ##store results
    system(paste("jellyfish dump", paste(qc_output, paste(basename(fasta), "counts.jf", sep="_"), sep="/"),">", paste(qc_output, paste(basename(fasta), "kmer_counts_dump.txt", sep="_"), sep="/"), sep=" "))
    ##graph kmer coverage
    ino.mers <- read.table(paste(qc_output, paste(basename(fasta), "kmer_hist.txt", sep="_"), sep="/"))
    jpeg(paste(output_dir, paste(basename(fasta), "kmer_histogram.jpg", sep="_"), sep="/"))
    kmer_histogram <- barplot(ino.mers[,2], ylab="No of kmers", xlab="Counts of a k-mer",
                              names.arg=ino.mers[,1], cex.names=0.8,
                              main=paste(basename(fasta), "Kmer Coverage", sep=" "))    
    dev.off()
    ## calc length distribution
    cat("Graphing length distribution...\n\n")
    ##read in fasta file
    my.fa <- readDNAStringSet(fasta)
    ##graph results
    par(las=2)
    jpeg(paste(output_dir, paste(basename(fasta), "length_dist.jpg", sep="_"), sep="/"))
    length_dist_jpg <-hist(width(my.fa), main=paste(basename(fasta), "Length Distribution", sep=" "),
                           ylab="Frequency",xlab="Length", breaks=100)
    dev.off()
}
##apply function to fasta files
lapply(fasta_files, run_QC)

##########################
## Blast hit coverage  ###
##########################
cat("Calculating blast hit coverage...\n\n")

##original assembly
system(paste(blast_cov, blastx_report, final.fa, blastdb.fa, paste(qc_output, paste(basename(final.fa), "blast_cov", sep="_"), sep="/"), ">", paste(qc_output, paste(basename(final.fa), "blast_cov", sep="_"), sep="/")))

##filtered assembly
system(paste(blast_cov, filtered_blastx, filtered.fa, blastdb.fa, paste(qc_output, paste(basename(filtered.fa), "blast_cov", sep="_"), sep="/"), ">", paste(qc_output, paste(basename(filtered.fa), "blast_cov", sep="_"), sep="/")))

##graph hit coverage
jpeg(paste(output_dir, "blast_hit_cov.jpg", sep="/"))
par(mfrow=c(1,2))
##read assembly coverage report
ass.cov <- read.table(file=paste(qc_output, paste(paste(basename(final.fa), "blast_cov", sep="_"), "w_pct_hit_length", sep="."), sep="/"), sep="\t")
##name coverage columns
colnames(ass.cov) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "db_hit_len", "pct_hit_len_aligned", "hit_descr")
##read filtered coverage report
filt.cov <- read.table(file=paste(qc_output, paste(paste(basename(filtered.fa), "blast_cov", sep="_"), "w_pct_hit_length", sep="."), sep="/"), sep="\t")
##name coverage columns
colnames(filt.cov) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "db_hit_len", "pct_hit_len_aligned", "hit_descr")
##graph assembly distribution
boxplot(ass.cov$pct_hit_len_aligned, main="Assembly blast hit coverage", ylab="% Coverage", xlab="Distribution", col="light green")
##graph filtered distribution
boxplot(filt.cov$pct_hit_len_aligned, main="Filtered blast hit coverage", ylab="% Coverage", xlab="Distribution", col="light blue")
dev.off()

##############################
## Bowtie Alignment Stats  ###
##############################
##run stats on bam file
system(paste("samtools flagstat", bam_file, ">", paste(qc_output, "bwt_align_stats.txt", sep="/"), sep=" "))

cat("\n Finished basic QC...\n\n")



 
