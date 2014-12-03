#!/usr/bin/env Rscript

###############################  trinnotato.R  ###########################################
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
dir.create(file.path(paste(output_dir, "trinotate", sep="/")))
trinnotate_out <- paste(output_dir, "trinotate", sep="/")

##change working directory to pep_output
setwd(trinnotate_out)

#######################
## Running trinotate ##
#######################
##function to run trinotate on fasta files
run_trinotate <- function(fasta){
  ##run hmmer
  system(paste(hmmer, "--cpu", CPU, "--domtblout", paste(basename(fasta), "Pfam.out", sep="_"),  pfam_path, paste(output_dir, "peptide", paste(basename(fasta), "transdecoder.pep", sep="."), sep="/"), 
               ">", paste(trinnotate_out, paste(basename(fasta), "pfam.out", sep="_"), sep="/")))
  ##run signalP
  system(paste(signalp, paste(output_dir, "peptide", paste(basename(fasta), "transdecoder.pep", sep="."), sep="/"), 
               ">", paste(trinnotate_out, paste(basename(fasta), "signalp.out", sep="_"), sep="/")))
  ##run tmHMM
  system(paste(tmHMM, paste(output_dir, "peptide", paste(basename(fasta), "transdecoder.pep", sep="."), sep="/"), 
               ">", paste(trinnotate_out, paste(basename(fasta), "tmHMM.out", sep="_"), sep="/")))
  ##make directory and setup sqlite database
  system(paste("mkdir", paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/")), wait=TRUE)
  system(paste("wget http://sourceforge.net/projects/trinotate/files/TRINOTATE_RESOURCES/20140708/Trinotate.20140708.swissprot.sqlite.gz/download -O", paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite.gz", sep="/")), wait=TRUE)
  system(paste("gunzip", paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite.gz", sep="/")), wait=TRUE)
  ##create gene to transcripts mapping file from RSEM output
  fasta_file <- readDNAStringSet(final.fa)
  ##create data frame of transcript id's
  fasta_names.df <- data.frame(transcript_id = unlist(mclapply(strsplit(as.character(names(fasta_file)), " "), function(x) x[1], mc.cores=5)))
  ##create data frame of gene id's
  fasta_names.df$gene_id <- unlist(mclapply(strsplit(as.character(fasta_names.df$transcript_id), "_"), function(x) paste(x[1], x[2], sep="_"), mc.cores=5))
  ##clean up scaffold gene id's
  fasta_names.df$gene_id <- gsub("_NA", "", fasta_names.df$gene_id)
  ##write to a file in correct order
  fasta_names.df <- fasta_names.df[,c(2,1)]
  write.table(fasta_names.df, file=paste(trinnotate_out, "genes_to_trans.txt", sep="/"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  ##load results into trinoate database
  setwd(trinnotate_output)
  ##load genes to trans mapping file
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "init --gene_trans_map", paste(trinnotate_out, "genes_to_trans.txt", sep="/"), "--transcript_fasta", final.fa, "--transdecoder_pep", paste(paste(output_dir, "peptide", paste(basename(fasta), "transdecoder.pep", sep="."), sep="/"))))  
  ##load blastp report
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "LOAD_swissprot_blastp", blastp_report))
  ##load blastx report
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "LOAD_swissprot_blastx", blastx_report))
  ##load Pfam report
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "LOAD_pfam", paste(trinnotate_out, paste(basename(fasta), "pfam.out", sep="_"), sep="/")))
  ##load tmHMM report
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "LOAD_tmhmm", paste(trinnotate_out, paste(basename(fasta), "tmHMM.out", sep="_"), sep="/")))
  ##load singalP predictions
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "LOAD_signalp", paste(trinnotate_out, paste(basename(fasta), "signalp.out", sep="_"), sep="/")))
  ##create annotation report
  system(paste(trinotate_path, paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), "report >", paste(trinnotate_out, paste(basename(fasta), "annotation.xls", sep="_"), sep="/")))
  ##load annotation report to database
  system(paste("/home/sre/tools/Trinotate_r20140708/util/annotation_importer/import_transcript_names.pl", paste(paste(trinnotate_out, paste(basename(fasta), "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), paste(trinnotate_out, paste(basename(fasta), "annotation.xls", sep="_"), sep="/")))
  ##copy results to project folder for viewing
  system(paste("mkdir", paste("/n/projects/sre/", paste(species, "sqlite_db", sep="_"), sep="/")), wait=TRUE)
  ##copy database to projects folder
  system(paste("cp", paste(paste(trinnotate_out, paste(species, "sqlite_db", sep="_"), sep="/"), "Trinotate.sqlite", sep="/"), paste("/n/projects/sre/", paste(basename(fasta), "sqlite_db/Trinotate.sqlite", sep="_"), sep="/")))
  ##copy needed trinotate folders into database directory
  system(paste("cp -Rf", paste(output_dir, "js", sep="/"), paste("/n/projects/sre/", paste(species, "sqlite_db/", sep="_"), sep="")))
  system(paste("cp -Rf", paste(output_dir, "PerlLib", sep="/"), paste("/n/projects/sre/", paste(species, "sqlite_db/", sep="_"), sep="")))
  system(paste("cp -Rf", paste(output_dir, "cgi-bin", sep="/"), paste("/n/projects/sre/", paste(species, "sqlite_db/", sep="_"), sep="")))
}

##apply function to final.fa
run_trinotate(final.fa)



