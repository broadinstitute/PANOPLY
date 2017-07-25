

source ('preamble.r')

bsub.file <- 'correlation.bsub'


system.x <- function (cmd, ...) {
  # deal with idiosyncracies of UGER at Broad
  # the environment is not inherited -- need to re-run 'use UGER' before running qsub/qstat
  pre.cmd <- "source /broad/software/scripts/useuse; reuse -q UGER"
  full.cmd <- sprintf ("%s; %s", pre.cmd, cmd)
  system (full.cmd, ...)
}


write.bsub.file <- function (prefix, pome) {
  script <- paste (prefix, '-', bsub.file, sep='')
  cat ("#!/bin/bash\n",
       "#$ -o ", prefix, "-dump/corrcalc-out-$TASK_ID.txt    # std out\n",
       "#$ -cwd                                 # run in current directory\n",
       "#$ -V                                   # propagate enrivonment\n",
       "#$ -j y                                 # combine stderr with stdout\n",
       "#$ -q long                              # queue\n",
       "#$ -tc 200                              # total concurrent jobs\n",
       "#$ -t 1-", LSF.mut.jid.max, "           # job array\n",
       "#$ -N ", prefix, "corr", "              # job array name\n\n",
       "source /broad/software/scripts/useuse \n",
       "reuse -q R-3.1 \n",
       file=script, sep='')
  cat ("Rscript correlation.r $SGE_TASK_ID", LSF.mut.jid.max, prefix, pome, "\n",
       file=script, append=TRUE)  
}


run.bsub.file <- function (prefix) {
  # remove dump and output directories and create new ones
  dump.dir <- paste (prefix, '-dump', sep='')
  out.dir <- paste (prefix, '-output', sep='')
  if (file.exists (dump.dir)) unlink (dump.dir, recursive=TRUE)
  if (file.exists (out.dir)) unlink (out.dir, recursive=TRUE)
  dir.create (dump.dir)
  dir.create (out.dir)
  
  cmd <- paste ('qsub < ', prefix, '-', bsub.file, sep='')
  system.x (cmd)
}


finish.corr <- function (prefix) { 
  # wait for all jobs to finish and then run the plotting program
  jobs.done <-FALSE
  while (! jobs.done) {
    jname <- paste (prefix, 'corr', sep='')
    jobs.done <- as.numeric (system.x (paste ('qstat -j', jname, '| wc -l'), intern=TRUE)) <= 1
    Sys.sleep (120)
  }
  
  cmd <- paste ('qsub -q long -o uger.out -cwd -j y -l h_vmem=16g -V runR.sh Rscript correlation.r 0 ', 
                length (list.files (paste (prefix, '-output', sep='')))/2, ' ', prefix, ' NULL', sep='')
  system.x (cmd)
}



prefix <- paste (type, subset.str, sep='')
write.bsub.file (prefix, pome=file.path (netgestalt.dir, paste (type, '-matrix.csv', sep='')))
run.bsub.file (prefix)
finish.corr (prefix)

