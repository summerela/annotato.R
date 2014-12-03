##file: program_config.R

## All variables must be enclosed in quotes ##

#############################################
############# program paths  ################
#############################################
## Trimming ##
trimmomatic <- "/home/sre/tools/Trinity07172014/trinity-plugins/Trimmomatic/trimmomatic.jar"

## Assembly ##
trinity <- "/home/sre/tools/Trinity07172014/Trinity"
SOAP <- "/n/local/bin/SOAPdenovo-Trans-31mer"

## merge fasta files ##
cd.hit.est <- "/home/sre/tools/cd-hit-v4.6.1-2012-08-27/cd-hit-est"
seqclean <- "/home/sre/tools/seqclean-x86_64/seqclean"

## Generate Base Files ##
pfam_path <- "/home/sre/tools/hmmer/pfam/Pfam-AB.hmm.bin"
rsem_path <- "/home/sre/tools/Trinity07172014/util/align_and_estimate_abundance.pl"

## Filtering programs ##
transPS <- "/home/sre/tools/TransPS1.1.0/transps.pl"
countFPKM <- "/home/sre/tools/Trinity07172014/util/misc/count_features_given_MIN_FPKM_threshold.pl"
##path to detonate distribution estimation script
detonate_dist <- "/home/sre/tools/detonate-1.8/rsem-eval/rsem-eval-estimate-transcript-length-distribution"
##path to detonate RSEM-eval score script
detonate <- "/home/sre/tools/detonate-1.8/rsem-eval/rsem-eval-calculate-score" 

## QC programs ##
##path to summary stats script
fastatools <- "/home/ejr/tools/ejrtools/fastatools"
blast_cov <- "/home/sre/tools/Trinity07172014/util/analyze_blastPlus_topHit_coverage.pl"
bwt_aln <- "/home/sre/tools/Trinity07172014/util/SAM_nameSorted_to_uniq_count_stats.pl"

## peptide programs ##
trans_path <- "/home/sre/tools/Trinity07172014/trinity-plugins/transdecoder/TransDecoder"

## trinnotate programs ##
hmmer <- "/home/sre/tools/hmmer/binaries/hmmscan"
signalp <- "/home/sre/tools/signalp-4.1/signalp"
tmHMM <- "/home/sre/tools/tmhmm-2.0c/bin/tmhmm"
RNAMMER <- "/home/sre/tools/Trinotate_r20140708/util/rnammer_support/RnammerTranscriptome.pl"
rnammer_path <- "/home/sre/tools/rnammer-1.2/rnammer"
trans_map <- "/home/sre/tools/Trinity07172014/util/support_scripts/get_Trinity_gene_to_trans_map.pl"
trinotate_path <- "/home/sre/tools/Trinotate_r20140708/Trinotate"