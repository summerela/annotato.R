#!/usr/bin/env Rscript

##########################################
##Setting up script enviorment ##
##########################################
##set global environment variables to locate packages
Sys.setenv(PATH="home/sre/.linuxbrew/bin:/n/local/stage/rbenv/rbenv-0.3.0/shims:/n/local/stage/rbenv/rbenv-0.3.0/shims:/n/local/stage/rbenv/rbenv-0.3.0/bin:/n/local/stage/quake/current/bin:/n/local/stage/pythonbrew/pythonbrew-1.0/bin:/n/local/stage/pythonbrew/pythonbrew-1.0/pythons/Python-2.7.5/bin:/n/local/stage/perlbrew/perlbrew-0.43/bin:/n/local/stage/perlbrew/perlbrew-0.43/perls/perl-5.16.1t/bin:/n/local/bin/mpich2:/n/local/bin/mpich2:/n/local/bin/meme:/n/local/bin:/n/local/bin/bio.brew/bin:/n/local/bin/bio.brew/bin:/n/local/bin:/n/local/proteomics/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/n/local/stage/go/current/bin:/n/local/stage/go/bin:/n/local/stage/go/bin:/n/local/stage/go/bin")

#each of the following scripts were run
scripts <- c("base_make.R", "filte.R", "stat_generato.R", "annotato.R", "translato.R")

##function to run scripts and output errors
assemble <- function(script) {
  out <- tryCatch(
{  
  readLines(con=scripts, warn=FALSE) 
},
error=function(cond) {
  message(paste("There was in error in", script, sep=" "))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  stop()
},
warning=function(cond) {
  message(paste("The", script, "script generated a warning:", sep=" "))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  stop()
},
finally={
  ##read in scripts
  source("user_config.R")
  source("script_config.R")
  source("program_config.R")
  ##run eaach script
  setwd(output_dir)
  cat(paste("\n Running", script, "\n", sep=" "))
  system(paste(output_dir, script, sep="/"), wait=TRUE)
  message(paste("Finished running:", script))
}
  )    
return(out)
}

##apply the function to each script in the list
lapply(scripts, assemble)


