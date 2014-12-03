#!/usr/bin/env Rscript
#############################################################################################
############################### assemble.R  ################################################# #############################################################################################
##  PURPOSE: This program takes as input transcript fastq files and assembles them
##  using Trinity at 25mer SOAPdenovo-trans at 19 and 31 kmers. The output fasta files
##  can be concatenated together and run through downstream analysis. It is recommended to
##  run CD-HIT-EST to remove redundancy between data sets.  
## 
##  Input: pre-trimmed (use trim_fastq.R) paired or single-end fastq files. Fastq files are
##         assumed to be located in the output directory specified when using trim_fastq.R.
##
##  Output:
##       Trinity.fa = trinity fasta file
##       SOAP19.fa = SOAPdenovo-trans scaffolded contigs at 19mers
##       SOAP31.fa.scafSeq = SOAPdenovo-trans scaffolded contigs at 31mers
##
##       The final output files will be placed in the output directory /assembly folder. 
################################################################################################
options(error=traceback)

##########################
## Read in config file  ##
##########################
source("user_config.R")
source("script_config.R")
source("program_config.R")

cat("\n Setting up output directory. \n")

##check that output directory exists
dir.create(file.path(output_dir))
##change working directory to output directory
setwd(output_dir)

##create directory for assembly
dir.create(file.path(paste(output_dir, "fasta_files", sep="/")))
assembly_out <- paste(output_dir, "fasta_files", sep="/")

##################################
##  Generate SOAP config file   ##
##################################
cat("Creating configuration file for assembly...\n")

##paired-end config file
if (reads == "paired"){
    line1 <- paste("max_rd_len=", max_read_length, sep="")
    line2 <- "[LIB]"
    line3 <- paste("rd_len_cutof=", as.integer(max_read_length)-5, sep="")
    line4 <- paste("avg_ins=", avg_ins, sep="")
    line5 <- paste("reverse_seq=", direction, sep="")
    line6 <- "asm_flags=3" ##align and scaffold
    line7 <- "map_len=32" ##keeping at default per manual
    line8 <- paste("q1=", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"), sep="") ##trimmed fastq location
    line9 <- paste("q2=", paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"), sep="")##trimmed fastq location
    ##merge together the commands, one per line
    config.file <- paste(line1, line2, line3, line4, line5, line6, line7, line8, line9, sep="\n")
    ##write the configuration file
    writeLines(config.file, paste(assembly_out, "paired_assembly.config", sep="/"))
}

##single-end config file
if (reads == "single") {
    line1 <- paste("max_rd_len=", max_read_length, sep="")
    line2 <- "[LIB]"
    line3 <- paste("rd_len_cutof=", as.integer(max_read_length)-5, sep="")
    line4 <- paste("avg_ins=", avg_read_length, sep="")
    line5 <- paste("reverse_seq=", direction, sep="")
    line6 <- "asm_flags=3"
    line7 <- "map_len=32"
    line8 <- paste("q=", paste(output_dir, "/trimmed/trimmed_single.fq", sep="/"), sep="")##trimmed fastq location
    config.file <- paste(line1, line2, line3, line4, line5, line6, line7, line8, sep="\n")
    writeLines(config.file, paste(assembly_out, "single_assembly.config", sep="/"))
}

##################################
##  Run only SOAPdenovo-trans   ##
##################################
##SOAP paired-end 
if (assembler == "2" & (reads == "paired"| reads == "PAIRED")){
    ##run the command
    SOAP19 <- paste(SOAP, "all -s", paste(assembly_out, "paired_assembly.config", sep="/"), "-o",
                    paste(assembly_out, "19mer.fa", sep="/"), "-K 19 -f -p", CPU, "-d 2", sep=" ")
    ##run SOAP 31mer
    SOAP31 <- paste(SOAP, "all -s", paste(assembly_out, "paired_assembly.config", sep="/"), "-o",
                 paste(assembly_out, "31mer.fa", sep="/"), "-K 31 -f -p", CPU, "-d 2", sep=" ")
}

##SOAP single-end 
if (assembler == "2" & reads == "single"| reads == "SINGLE"){
    ##run the command
    SOAP19 <- paste(SOAP, "all -s", paste(output_dir, "single_assembly.config", sep="/"), "-o",
                    paste(assembly_out, "19mer.fa", sep="/"), "-K 19 -f -p", CPU, "-d 2", sep=" ")
    SOAP31 <- paste(SOAP, "all -s", paste(output_dir, "single_assembly.config", sep="/"), "-o",
                  paste(assembly_out, "19mer.fa", sep="/"), "-K 31 -f -p", CPU, "-d 2", sep=" ")
}

#############################
##  Run only Trinity 25mer ##
#############################
##paired-end
if (assembler == "1" & (reads == "paired"| reads=="PAIRED")){
    ##create directionality options for trinity from user args 
    if (direction=="1"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --left", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"),
                           "--right",  paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"), "--SS_lib_type RF --CPU",
                           CPU, "--full_cleanup --min_kmer_cov 2 --bflyCalculateCPU --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="0"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --left", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"),
                           "--right",  paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"),
                           "--SS_lib_type FR --CPU", CPU, "--full_cleanup --min_kmer_cov 2 --bflyCalculateCPU --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="2"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --left", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"),
                           "--right",  paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"), "--CPU", CPU,
                           "--full_cleanup --min_kmer_cov 2 --bflyCalculateCPU --inchworm_cpu", CPU, sep=" ")
    }
}
##single-end
if (assembler == "1" & (reads == "single"| reads=="SINGLE")){
    ##create directionality options for trinity from user args 
    if (direction=="1"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --single",
                           paste(trimmed_output,"trimmed_single.fq", sep="/"), "--SS_lib_type RF --CPU", CPU,
                           "--full_cleanup --min_kmer_cov 2 --bflyCalculateCPU --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="0"){
        ##run trinity command
      Trinity25 <- paste(trinity, "--seqType fq --JM 100G --single",
                           paste(trimmed_output,"trimmed_single.fq", sep="/"), "--SS_lib_type FR --CPU", CPU,
                           "--full_cleanup --min_kmer_cov 2 --bflyCalculateCPU --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="2"){
        ##run trinity command
       Trinity25 <- paste(trinity, "--seqType fq --JM 100G --single",
                           paste(trimmed_output,"trimmed_single.fq", sep="/"), "--CPU", CPU,
                          "--full_cleanup --min_kmer_cov 2 --bflyCalculateCPU --inchworm_cpu", CPU, sep=" ")
   }
}

###########################
## Run SOAP and Trinity ##
##########################
##run both assemblers on paired end reads
if (assembler == "3" & (reads == "paired"| reads == "PAIRED")){
    ##run SOAP 19mer
    SOAP19 <- paste(SOAP, "all -s", paste(assembly_out, "paired_assembly.config", sep="/"), "-o",
                 paste(assembly_out, "19mer.fa", sep="/"), " -K 19 -f -p", CPU, "-d 2", sep=" ")

    ##run SOAP 31mer
    SOAP31 <- paste(SOAP, "all -s", paste(assembly_out, "paired_assembly.config", sep="/"), "-o",
                 paste(assembly_out, "31mer.fa", sep="/"), "-K 31 -f -p", CPU, "-d 2", sep=" ")

    ##create directionality options for trinity from user args 
    if (direction=="1"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --left", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"),
                           "--right",  paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"),
                           "--SS_lib_type RF --CPU", CPU,
                           "--full_cleanup --min_kmer_cov 2 --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="0"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --left", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"),
                           "--right",  paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"),
                     "--SS_lib_type FR --CPU", CPU,
                     "--full_cleanup --min_kmer_cov 2 --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="2"){
        ##run trinity command
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --left", paste(output_dir, "/trimmed/trimmed_left.fq", sep="/"),
                           "--right",  paste(output_dir, "/trimmed/trimmed_right.fq", sep="/"), "--CPU", CPU,
                           "--full_cleanup --min_kmer_cov 2 --inchworm_cpu", CPU, sep=" ")
    }
}


##Run both assemblers on single-end reads
if (assembler == "3" & (reads == "single"| reads == "SINGLE")){
    ##run the 19mer command
    SOAP19 <- paste(SOAP, "all -s", paste(assembly_out, "single_assembly.config", sep="/"), "-o",
                 paste(assembly_out, "19mer.fa", sep="/"), " -K 19 -f -p", CPU, "-d 2", sep=" ")

    ##run the 31mer command
    SOAP31 <- paste(SOAP, "all -s", paste(assembly_out, "single_assembly.config", sep="/"), "-o",
                  paste(assembly_out, "31mer.fa", sep="/"), "-K 31 -f -p", CPU, "-d 2", sep=" ")
    ##determine directionality
    if (direction=="1"){
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --single", paste(output_dir, "/trimmed/trimmed_single.fq", sep="/"),
                           "--SS_lib_type RF --CPU", CPU, "--full_cleanup --min_kmer_cov 2", sep=" ")
    }
    if (direction=="0"){
        Trinity25 <- paste(trinity, "--seqType fq --JM 100G --single",
                           paste(output_dir, "/trimmed/trimmed_single.fq", sep="/"), "--SS_lib_type FR --CPU",
                 CPU, "--full_cleanup --min_kmer_cov 2 --inchworm_cpu", CPU, sep=" ")
    }
    if (direction=="2"){
        Trinity <- paste(trinity, "--seqType fq --JM 100G --single", paste(output_dir, "/trimmed/trimmed_single.fq", sep="/"),
                         "--CPU", CPU, "--full_cleanup --min_kmer_cov 2 --inchworm_cpu", CPU, sep=" ")
    }
}

#########################
## create qsub command ##
#########################
qsub_cmd <- paste("qsub -N Assemble -j y -pe by_node_32", CPU, "-p -100 -l mem_free=4G,h_vmem=4.2G -M", email, "-m abe -b y", sep=" ")

#####################
## submit command ##
#####################
if (qsub == "true" | qsub == "TRUE"){
    cat("\n Running SOAP 19mer...\n")
    system(paste(qsub_cmd, SOAP19, sep=" "))
    cat("\n Running SOAP 31mer...\n")
    system(paste(qsub_cmd, SOAP31, sep=" "))
    cat("\n Running Trinity...\n")
    system(paste(qsub_cmd, Trinity25, sep=" "))
}
if (qsub == "false" | qsub == "FALSE"){
    cat("\n Running SOAP 19mer...\n")
    system(paste(SOAP19, sep=" "))
    cat("\n Running SOAP 31mer...\n")
    system(paste(SOAP31, sep=" "))
    cat("\n Running Trinity...\n")
    system(paste(Trinity25, sep=" "))
}

###################
## run cleanup  ###
###################
cat("\n Cleaning up.. \n")
## rename fasta files
system(paste("mv", paste(assembly_out, "19mer.fa.scafSeq", sep="/"), paste(assembly_out, "SOAP19.fa", sep="/"), sep=" "))
system(paste("mv", paste(assembly_out, "31mer.fa.scafSeq", sep="/"), paste(assembly_out, "SOAP31.fa", sep="/"), sep=" "))
system(paste("mv", paste("*.Trinity.fasta", paste(assembly_out, "Trinity.fa", sep="/")), sep=" "))

##remove 19 and 31 mer messsy files
system(paste("rm -Rf", paste(assembly_out, "19mer.*", sep="/"), sep=" "))
system(paste("rm -Rf", paste(assembly_out, "31mer.*", sep="/"), sep=" "))

       
cat("\nAssembly(s) completed. Please check stdout for any error messages.\n")
