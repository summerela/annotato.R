##user_config.R

#############################################################################
## update the paths for each of the following variables. All lists         ##
## should be enclosed in quotes, and separated by a comma with no spaces.  ## 
#############################################################################
##what is the name of the species being assembled (for report)
species <- "Test Lizard"

##directory path(s) to fastq files
fastq_dir <- "/home/sre/testing/"

##specify if your reads are paired or single
reads <- "paired"

##Enter 3 for TruSeq3 (as used by HiSeq and MiSeq machines)
##or 2 for TruSeq2 (as used in GAII machines)
illumina <- "3"

##enter the adapters used for creating your libraries
adapters <- "ACAGTG,AGGATT"

##enter the max read length for your reads
##get info from mol bio
max_read_length = "151"

##enter the average insert size of your reads
##get info from mol bio
avg_ins = "450"

##enter stranded direction
## 0 = FR, 1 = RF (dUTP), 2 = non-directional
direction = "1"

##select which assemblers to run
## 1 = Trinity only, 2 = SOAP only, 3 = both
assembler = "3"

##specify number of cores to use for each job (cluster can only be 16 or 32)
CPU <- "12"

##Specify "true" if submitting the job to the cluster, false for a regular server
qsub <- "false"

##specify an email address for the cluster to send job notiications
email <- "sre@stowers.org"

##specify the full file path where you want to place the output files
output_dir <- "/home/sre/studio_test/"

##Contigs with at least this % nucleotide identity will be merged
##recommend starting with 0.99
PID <- "0.99"

##SeqClean memory limit (in MB); 0 for unlimited (recommend 60000)
memory <- "60000"

##enter location of protein database used for blast
##include everything before the .p
blastdb <- "/home/sre/BlastDB/uni_sp_prot"

##select isoform or gene to specify if you want to collapse results by gene
##or keep all isoforms that meet filtering criteria
filter_method <- "isoform"

##path to fasta file  used to create blast database, such as the swissprot.fa
blastdb.fa <- "/home/sre/BlastDB/fasta/uniprot_sprot.fasta"

##enter full path to a nucelotide fasta of the closest related annotated species
ref.fa <- "/home/sre/BlastDB/fasta/Anolis_carolinensis.AnoCar2.0.75.cdna.all.fa"

##enter evalue cutoff for blast hits
eval <- ".01"

##enter minimum percent identity to blast hit required
blast_PID <- "30"

##enter minimum base pairs of match length
min_blast_length <- "50"

