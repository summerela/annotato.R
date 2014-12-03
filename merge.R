#!/usr/bin/env Rscript
###############################################################################
###############          merge.R    ########################################### ###############################################################################
##  PURPOSE: This program takes as input fasta files from different assemblies
##           and removes redundancy from each file by running CD-HIT-EST
##           to filter all contigs that match with at least 99% (or other
##           user-specified amount) identity, keeping the longer contig. The
##           program then filters out transcripts shorter than 200 bp,removes
##           transcripts that are more than 99% identical on a nucelotide level
##           to human, staph, ecoli, pombe, cerevisiae,drosophila and univec
##           sequences. The assemblies are then merged into a final,
##           non-redundant fasta.
##
##  Input: takes a list of fasta files to merge 
##
##  Output:
##       fasta_200 = fasta after removing transcripts < 200 bp
##       fasta_reduced = fasta_200 after removing redundancy
##       merged.fa = fasta file containing the merged reads from each fasta file
##                   after removing <200 bp and redundant reads 
##       final.fa = final fasta after merge and clean
##
## NOTE: This program will merge all fasta files ending in "_reduced" in the working
##       directory! You should also update the program and vector/contaminant
##       database paths for your system. 
###################################################################################
##check that packages are installed, and if not, install
packages <- c("multicore", "Biostrings")
##install any packages that are not currently installed
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.r-project.org")
}
##load packages and supress their startup messages
suppressPackageStartupMessages(lapply(packages, require, character.only=T))

##########################
## Read in config file  ##
##########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

##check that output directory exists
dir.create(file.path(output_dir))

##check if output directory exists, if not, create it
dir.create(file.path(paste(output_dir, "merged", sep="/")))
merged_output <- paste(output_dir, "merged", sep="/")

##change working directory to output directory
setwd(merged_output)

####################################
## Gather and filter fasta files  ##
####################################
## Read in fasta list ##
fasta.list <- list.files(path=paste(output_dir, "fasta_files", sep="/"), pattern="*.fa|*.fasta", full.names=TRUE, ignore.case=TRUE)

##Filter out transcripts < 200 bp
min200 <- function(x){
    f <- readDNAStringSet(x)
    f1 <- f[which(width(f)>200)]
    writeXStringSet(f1, file=paste(paste(merged_output, basename(x), sep="/"), "200", sep="_"))
}
##apply function to each file in fasta.list
lapply(fasta.list, min200)

##create list of filtered files
fasta200 <- do.call("rbind", lapply(merged_output, list.files, pattern="*_200"))

cd.hit <- function(fa){
    ##create CD-HIT-EST command
    cd.hit_cmd <- paste(cd.hit.est, "-T", CPU, "-i", paste(merged_output, fa, sep="/"), "-c", PID, "-o",
                        paste(merged_output, paste(basename(fa), "reduced", sep="_"), sep="/"),
                        "-M", memory, sep=" ")
     system(cd.hit_cmd, wait=TRUE)
}

##############################
## create SeqClean command  ##
##############################
clean <- paste(seqclean, paste(merged_output, "final_merged.fa", sep="/"), "-o", paste(merged_output, "final_merged.fa.clean", sep="/"), "-v /home/sre/BlastDB/fasta/UniVec.fa -l 200 -c 16 -x 0.99 -s /home/ejr/db/human_genome.nt,/home/ejr/db/Dmel_5_48.nt,/home/ejr/db/ecoli.nt,/home/ejr/db/bacterial_16s.nt,/home/sre/BlastDB/staph.fa,/home/sre/BlastDB/cerevisiae.fa,/home/sre/BlastDB/pombe.fa", sep=" ")

#########################
## create qsub command ##
#########################
qsub_cmd <- paste("qsub -N Merge_fa -j y -pe by_node_32", CPU, "-p -100 -l mem_free=4G,h_vmem=4.2G -M", email, "-m abe -b y", sep=" ")

#####################
## submit command ##
#####################
##submit to cluster
if (qsub == "true" | qsub == "TRUE"){
    cat("\n \n Removing redundancy from fasta ...\n \n")
    ##prevent cd-hit memory errors
    system(paste("ulimit -s unlimited"))
    ##run CD-HIT-EST on each fasta
    system(paste(qsub_cmd, mclapply(fasta200, cd.hit, mc.cores=5), sep=" "), wait=TRUE)
    ## merge files
    cat("\n \n Merging fasta files...\n \n")
    system(paste(qsub_cmd, "cat", paste(paste(merged_output, list.files(merged_output, "_reduced$"), sep="/"), collapse=" "), ">", paste(merged_output, "merged.fa", sep="/"), sep=" "), wait=TRUE)
    ##run cd-hit on merged fasta to make seqclean faster
    cat("\n \n Removing redundancy...\n \n")
    system(paste(qsub_cmd, cd.hit.est, "-T", CPU, "-i", paste(merged_output, "merged.fa", sep="/"), "-c", PID,
                 "-o", paste(merged_output, "final_merged.fa", sep="/"), "-M", memory, sep=" "), wait=TRUE)
    ##run seqclean on merged files
    cat("\n \n Running SeqClean... \n \n")
    system(paste(qsub_cmd, clean, sep=" "), wait=TRUE)
    ##remove cleaning folders
    cat("\n \n Cleaning up... \n \n")
    system(paste("rm -Rf cleaning_*"))
    ##run cd-hit on cleaned fasta because it somehow introduces redundancy
    cat("\n \n Merging final output...\n \n")
    system(paste(qsub_cmd, cd.hit.est, "-T", CPU, "-i", paste(merged_output, "final_merged.fa.clean", sep="/"), "-c", PID,
                 "-o", paste(merged_output, "final.fa", sep="/"), "-M", memory, sep=" "))
    cat("\n \n Final output will be written to output_dir/merged/final.fa...\n \n")
}
##submit to server
if (qsub == "false" | qsub == "FALSE"){
   cat("\n \n Removing redundancy from fasta ...\n \n")
    ##prevent cd-hit memory errors
    system(paste("ulimit -s unlimited"))
    ##run CD-HIT-EST on each fasta
    system(paste(mclapply(fasta200, cd.hit, mc.cores=5)), wait=TRUE)
    ## merge files
    cat("\n \n Merging fasta files...\n \n")
    system(paste("cat", paste(paste(merged_output, list.files(merged_output, "_reduced$"), sep="/"), collapse=" "), ">", paste(merged_output, "merged.fa", sep="/"), sep=" "), wait=TRUE)
    ##run cd-hit on merged fasta
    cat("\n \n Merging final output...\n \n")
    system(paste(cd.hit.est, "-T", CPU, "-i", paste(merged_output, "merged.fa", sep="/"), "-c", PID,
                 "-o", paste(merged_output, "final_merged.fa", sep="/"), "-M", memory, sep=" "), wait=TRUE)
   ##run seqclean on merged files
   cat("\n \n Running SeqClean...\n \n")
   system(paste(clean), wait=TRUE)
   ##remove cleaning folders
   cat("\n \n Cleaning up... \n \n")
   system(paste("rm -Rf cleaning_*"))
    ##run cd-hit on cleaned fasta because it somehow introduces redundancy
    cat("\n \n Merging final output...\n \n")
    system(paste(cd.hit.est, "-T", CPU, "-i", paste(merged_output, "final_merged.fa.clean", sep="/"), "-c", PID,
                 "-o", paste(merged_output, "final.fa", sep="/"), "-M", memory, sep=" "))
   cat("\n \n Complete. Final output written to ouptut_dir/merged/final.fa...\n \n")
}


