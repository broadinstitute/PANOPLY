

# check to make sure libraries needed by correlation.r are already installed
# (else, many jobs can fail before one of the jobs ends up installing the libraries)
if (!require("pacman")) install.packages ("pacman")
pacman::p_load (psych)


source ('preamble.r')

bsub.file <- 'correlation.bsub'


write.bsub.file <- function (prefix, pome) {
  script <- paste (prefix, '-', bsub.file, sep='')
  switch (compute.cluster.type,
          uger={
            # Sun Grid Engine (SGE/UGER)
            cat ("#!/bin/bash\n",
                 "#$ -o ", scratch.fs, prefix, "-dump/corrcalc-out-$TASK_ID.txt    # std out\n",
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
          },
          slurm={
            # SLURM (see https://srcc.stanford.edu/sge-slurm-conversion for SGE to SLURM conversion)
            cat ("#!/bin/bash\n",
                 "#SBATCH -o ", scratch.fs, prefix, "-dump/corrcalc-out-%a.txt    # std out\n",
                 "#SBATCH --export=ALL                                   # propagate enrivonment\n",
                 "#SBATCH -n 1                                           # tasks, 1 for R\n",
                 "#SBATCH -p ", cluster.queue, "                         # partition/queue to use\n",
                 "#SBATCH -t 8:0:0                                       # run time, required on some systems\n",
                 "#SBATCH --mem=8G                                       # memory \n",
                 "#SBATCH -a 1-", LSF.mut.jid.max, "                     # job array\n",
                 "#SBATCH -J ", prefix, "corr", "                        # job array name\n\n",
                 file=script, sep='')
            cat ("Rscript correlation.r $SLURM_ARRAY_TASK_ID", LSF.mut.jid.max, prefix, pome, "\n",
                 file=script, append=TRUE)
          },
          {
            # default -- run serially; may crash or fail
            cat ("#!/bin/bash\n", file=script, sep='')
            cat ("Rscript correlation.r 1 1", prefix, pome, "\n", file=script, append=TRUE)
          }
  )
}


run.bsub.file <- function (prefix) {
  # remove dump and output directories and create new ones
  dump.dir <- paste (scratch.fs, prefix, '-dump', sep='')
  out.dir <- paste (scratch.fs, prefix, '-output', sep='')
  if (file.exists (dump.dir)) unlink (dump.dir, recursive=TRUE)
  if (file.exists (out.dir)) unlink (out.dir, recursive=TRUE)
  mkdir (dump.dir)
  mkdir (out.dir)
  if (scratch.fs != "") {
    # create links from cwd -- will generate warning if links exist
    file.symlink (dump.dir, ".")
    file.symlink (out.dir, ".")
  }
  
  cmd <- switch (compute.cluster.type,
                 uger=paste ('qsub < ', prefix, '-', bsub.file, sep=''),
                 slurm=paste ('sbatch ', prefix, '-', bsub.file, sep=''),
                 paste ('sh ', prefix, '-', bsub.file, sep='')
  )
  system.x (cmd)
}


finish.corr <- function (prefix) { 
  # wait for all jobs to finish and then run the plotting program
  jobs.done <-FALSE
  while (! jobs.done) {
    jname <- paste (prefix, 'corr', sep='')
    jobs.done <- compute.job.done (jname)
    Sys.sleep (120)
  }
  
  cmd <- switch (compute.cluster.type,
                 uger={paste ('qsub -q long -o uger.out -cwd -j y -l h_vmem=16g -V runR-uger.sh Rscript correlation.r 0 ', 
                              length (list.files (paste (prefix, '-output', sep='')))/2, ' ', prefix, ' NULL', sep='')},
                 slurm={paste ('sbatch -t 2:0:0 -p ', cluster.queue, ' -o slurm.out --mem=16G --export=ALL runR-slurm.sh Rscript correlation.r 0 ', 
                               length (list.files (paste (prefix, '-output', sep='')))/2, ' ', prefix, ' NULL', sep='')},
                 {paste ('Rscript correlation.r 0 1 ', prefix, ' NULL', sep='')}
  )
  system.x (cmd)
}



prefix <- paste (type, subset.str, sep='')
write.bsub.file (prefix, pome=file.path (netgestalt.dir, paste (type, '-matrix.csv', sep='')))
run.bsub.file (prefix)
finish.corr (prefix)

