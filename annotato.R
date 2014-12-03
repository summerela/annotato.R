#!/usr/bin/env Rscript


###############################  annotato.R  #################################################
##  PURPOSE: This program takes as input full and filtered nucleotide fastas and corresponding
##           blast reports and bam files. The ouput is an annotated fasta file consisting
##           of original contig name, gene description from best blast hit, and RPKM.
##           The output is also saved as an r dataframe for downstream analysis. 
##  OUTPUT:
##        fasta.fa_counts = read abundance counts, can be used for DE analysis
##        fasta.fa_annot = a fasta file with annotated headers, in the format:
##                    contig_name gi_accession gene_descrip RPKM
##        fasta.fa_annot.df = an R data frame of the annotation for downstream analysis
##        fasta.fa_nucDB.n* = set of three blast nucleotide database files
##
##        The final output file will be placed in the output directory specified in user args,
##        and the original fasta file name will be suffixed with _annot. 
###############################################################################################
options(error=traceback)
########################################
## Setting up libraries and variables ##
########################################
##check that packages are installed, and if not, install
install.packages(c("Biostrings", "edgeR", "multicore", "CHNOSZ"))

#########################
## Setting up script   ##
##########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

##check that output directory exists
dir.create(file.path(output_dir))

##create base_files directory
dir.create(file.path(paste(output_dir, "annotated", sep="/")))
annot_output <- paste(output_dir, "annotated", sep="/")

##create annotation function
annotate <- function(fasta, blast){
    ##change working directory to output directory
    setwd(annot_output)
    cat("\n Creating counts table...\n")
     ## Create counts table  ##
    system(paste("samtools idxstats", bam_file, ">", paste(annot_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "counts", sep="_"), sep="/"), sep=" "))
    ##read in the counts table
    counts <- read.table(paste(annot_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "counts", sep="_"), sep="/"), sep="\t")
    colnames(counts) <- c("contig", "length", "map", "unmap")
    ## Calculate RPKM  ##
    cat("\n Calculating RPKM... \n")
    ##read in blast file, filtering by user args set in config file
    blast_report <- read.blast(file=blast, similarity=blast_PID, evalue=eval, max.hits=1, min.length=min_blast_length)
    ##RPKM calc'ed as= (10^9 * Number of reads mapped to gene)/(Total Mapped Reads * Length of transcript)
    counts$RPKM <- rpkm(counts$map, gene.length=counts$length)
    ##  Add Sequences  ##
    cat("\n Reading in the fasta file...\n")
    ##read in fasta file
    fasta.fa <- readDNAStringSet(fasta)
    ##remove anything after the first space in fasta headers to make annotation pretty
    names(fasta.fa) <- mclapply(strsplit(names(fasta.fa), " "), function(x) x[1], mc.cores=5)
    ##grab sequences 
    counts$seqs <- as.character(fasta.fa)[match(as.character(counts$contig),as.character(names(fasta.fa)))]
    ##remove unmapped reads that contain no sequences in fasta file
    counts <- counts[!is.na(counts$seqs),]
    ## Add gene descriptions and blast info to dataframe  ##
    ##load in blast database subject fasta
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
    seqs <- DNAStringSet(counts$seqs)
    ##add annotation
    names(seqs) <- counts$annot
    ##save results to a fasta file
    writeXStringSet(seqs, filepath=paste(annot_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "annot", sep="_"), sep="/"))
    ##save counts object as r data object for future use
    save(counts, file=paste(annot_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "annot.df", sep="_"), sep="/"))
    ## create blast database from results
    cat("\n Creating blast database...\n")
    system(paste("makeblastdb -in", paste(annot_output, paste(mclapply(strsplit(as.character(fasta), "\\/"), function(x) x[length(x)], mc.cores=5), "annot", sep="_"), sep="/"), "-input_type fasta -dbtype nucl -out", paste(annot_output, paste(basename(fasta), "nucDB", sep="_"), sep="/"), sep=" "))    
    ##notify user that annotation has finished
    cat(paste("\n Finished annotating fasta files.", sep=" "))
}

##create data frame of fasta and blast reports
fasta_files<- c(final.fa, filtered.fa)
blast_reports <- c(blastx_report, filtered_blastx)
reports.df <- data.frame(fasta=unlist(as.character(fasta_files)), blasts=unlist(as.character(blast_reports)))

##apply annotation function to fasta list
mapply(annotate, fasta=as.character(reports.df$fasta), blast=as.character(reports.df$blasts))
 


