#!/usr/bin/env Rscript

###############################  translato.R  ###########################################
##  PURPOSE: This program takes as input full and filtered nucleotide fasta and bam, and 
##           translates them into peptide sequences. The peptide fastas are then blasted
##           against the user specified blast database and annoated with RPKM and
##           gene description matching its best hit. 
##  OUTPUT:
##       peptide.fa = nucleotide fasta translated into protein sequences
##               You will also be provided with corresponding bed, gff3 and cds files
##       peptide.df = R data frame of annotation
##       fasta_blastp = blastp report of protein query against protein db, tab separated
##       pepDB.p* = peptide blast database files
#########################################################################################
options(error=traceback)
########################################
## Setting up libraries and variables ##
########################################
##check that packages are installed, and if not, install
install.packages(c("Biostrings", "edgeR", "multicore", "CHNOSZ"))

#########################
## Setting up script   ##
#########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

##check that output directory exists
dir.create(file.path(output_dir))

##create base_files directory
dir.create(file.path(paste(output_dir, "peptide", sep="/")))
pep_output <- paste(output_dir, "peptide", sep="/")

##change working directory to pep_output
setwd(pep_output)

## create qsub command ##
#########################
qsub_cmd <- paste("qsub -N Translate -j y -pe by_node_32", CPU, "-p -100 -l mem_free=4G,h_vmem=4.2G -M", email, "-m abe -b y", sep=" ")

##list fastas from config files
fasta_files<- c(final.fa, filtered.fa)

################################
## Run Transdecoder and Blast ##
###############################
translate <- function(fasta){
    ##create transdecoder command
    if (direction == "O"){
        trans_cmd <- paste(trans_path, "-t", fasta, "--CPU", CPU, "-S -m 50 --search_pfam", pfam_path, sep=" ")
    }
    if (direction == "1"){
        trans_cmd <- paste(trans_path, "-t", fasta, "--CPU", CPU, "-S -m 50 --search_pfam", pfam_path, sep=" ")
    }
    if (direction == "2"){
        trans_cmd <- paste(trans_path, "-t", fasta, "--CPU", CPU, "-m 50 --search_pfam", pfam_path, sep=" ")
    }
    ##submit command to cluster or server
    if (qsub == "true" | qsub == "TRUE"){
        cat("\n Submitting Transdecoder...\n")
        system(paste(qsub_cmd, trans_cmd, sep=" "), wait=TRUE)
    }
    if (qsub == "false" | qsub == "FALSE"){
        cat("\n Running Transdecoder...\n")
        system(paste(trans_cmd), wait=TRUE)
    }
    ##create blastp command
    cat("\n Running blastp...\n")
    blastp_cmd <- paste("blastp -db", blastdb, "-query", paste(pep_output, paste(basename(fasta),"transdecoder.pep", sep="."),sep="/"),"-outfmt 6 -num_threads", CPU, "-evalue .01 -max_target_seqs 1 -out",paste(pep_output, paste(basename(fasta), "blastp", sep="_"), sep="/"))
    ##submit command to cluster or server
    if (qsub == "true" | qsub == "TRUE"){
        cat("\n Submitting blastp...\n")
        system(paste(qsub_cmd, blastp_cmd, sep=" "), wait=TRUE)
    }
    if (qsub == "false" | qsub == "FALSE"){
        cat("\n Running blastp...\n")
        system(paste(blastp_cmd), wait=TRUE)
    }
}
##run transdecoder and blastp on fasta files
lapply(fasta_files, translate)

##read in blast reports
pep_files <- list.files(path=pep_output, pattern="*.transdecoder.pep$", full.names=TRUE)
blast_reports <- list.files(path=pep_output, pattern="*_blastp$", full.names=TRUE)
reports.df <- data.frame(fasta=unlist(as.character(pep_files)), blasts=unlist(as.character(blast_reports)))

##create annotation function
annotate <- function(fasta, blast){
    ##change working directory to output directory
    setwd(pep_output)
    cat("\n Creating counts table...\n")
     ## Create counts table  ##
    system(paste("samtools idxstats", bam_file, ">", paste(pep_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "counts", sep="_"), sep="/"), sep=" "))
    ##read in the counts table
    counts <- read.table(paste(pep_output, paste(mclapply(strsplit(as.character(fasta), "/"), function(x) x[length(x)], mc.cores=5), "counts", sep="_"), sep="/"), sep="\t")
    colnames(counts) <- c("contig", "length", "map", "unmap")    
    ##remove unmapped reads 
    counts <- counts[which(counts$map != 0),]
    ## Calculate RPKM  ##
    cat("\n Calculating RPKM... \n")
    ##read in blast file, filtering by user args set in config file
    blast_report <- read.blast(file=blast, similarity=blast_PID, evalue=eval, max.hits=1)
    ##edit blast qseqid to match with counts and fasta
    blast_report$queryId <- mclapply(strsplit(as.character(blast_report$queryId), "\\|"), function(x) x[1], mc.cores=5)
    ##RPKM calc'ed as= (10^9 * Number of reads mapped to gene)/(Total Mapped Reads * Length of transcript)
    counts$RPKM <- rpkm(counts$map, gene.length=counts$length, normalized.lib.sizes=TRUE)
    ##  Add Sequences  ##
    cat("\n Reading in the fasta file...\n")
    ##read in fasta file
    fasta.fa <- readAAStringSet(fasta)
    ##remove anything after the first space in fasta headers to make annotation pretty
    names(fasta.fa) <- mclapply(strsplit(names(fasta.fa), "\\|"), function(x) x[1], mc.cores=5)
    ##grab sequences 
    counts$seqs <- as.character(fasta.fa)[match(as.character(counts$contig),as.character(names(fasta.fa)))]
    ##remove unmapped reads 
    counts <- counts[which(counts$seqs != "NA"),]
    ## Add gene descriptions and blast info to dataframe  ##
    ##load in blast database subject fasta
    cat("\n Reading in blast database fasta...\n")
    subject.fa <- readAAStringSet(blastdb.fa)
    ##data frame subject fasta
    subject.df <- data.frame(comp=unlist(mclapply(strsplit(as.character(names(subject.fa)), " "), function(x) x[1], mc.cores=5)), descrip=unlist(mclapply(strsplit(as.character(names(subject.fa)), " "), function(x) paste(x[-1], collapse=" "), mc.cores=5)))
    cat("\n Merging blast, counts and sequence information...\n")
    ##add blast subjectID to counts dataframe
    counts$subject <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$subjectId    
    ##add blast PID counts to dataframe
    counts$PID <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$percIdentity
    ##add blast alnLength to counts dataframe
    counts$match_length <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$alnLength
    ##add blast  queryStart to counts dataframe
    counts$qstart <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$queryStart
    ##add blast  queryStart to counts dataframe
    counts$qend <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$queryEnd
    ##add blast  queryStart to counts dataframe
    counts$sstart <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$subjectStart
    ##add blast  queryStart to counts dataframe
    counts$send <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$subjectEnd
    ##add blast eVal to counts dataframe
    counts$eval <- blast_report[match(as.character(counts$contig), as.character(blast_report$queryId)),]$eVal
    ##add swissp gene descriptions to counts
    counts$descrip <- subject.df[match(counts$subject, subject.df$comp),]$descrip
    cat("\n Annotating your sequences...\n")
    ## Create Annotations  ##
    counts$annot <- paste(as.character(counts$contig), as.character(counts$descrip),
                          "RPKM", formatC(counts$RPKM, digits=2, format="f"), sep=" ")
    ##remove software breaking characters from annotation names
    counts$annot <- chartr("!*@#$%^:=", "__________", counts$annot)
    ## Save Results ##
    cat("\n Saving results... \n")
    ##create fasta file format
    seqs <- AAStringSet(counts$seqs)
    ##add annotation
    names(seqs) <- counts$annot
    ##save results to a fasta file
    writeXStringSet(seqs, filepath=paste(pep_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "annot", sep="_"), sep="/"))
    ##save counts object as r data object for future use
    save(counts, file=paste(pep_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "annot.df", sep="_"), sep="/"))
    ## create blast database from results
    cat("\n Creating blast database...\n")
    system(paste("makeblastdb -in", paste(pep_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "annot", sep="_"), sep="/"), "-input_type fasta -dbtype prot -out", paste(pep_output, paste(basename(fasta), "pepDB", sep="_"), sep="/"), sep=" "))    
    ##notify user that annotation has finished
    cat(paste("\n Finished annotating peptide file.\n", sep=" "))
}

##apply annotation function to fasta list
mapply(annotate, fasta=as.character(reports.df$fasta), blast=as.character(reports.df$blasts))

