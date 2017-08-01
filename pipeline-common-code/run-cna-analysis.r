
source ('preamble.r')
min.cna.N <- 5
bsub.file <- 'cna-analysis.bsub'


script.name <- function (prefix) paste (prefix, '-', bsub.file, sep='')


write.bsub.file <- function (prefix, mrna, cna, pome) {
  script <- script.name (prefix)
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
            cat ("Rscript cna-analysis.r $SGE_TASK_ID", LSF.mut.jid.max, prefix, mrna, cna, pome, "\n",
                 file=script, append=TRUE)  
          },
          slurm={
            # SLURM (see https://srcc.stanford.edu/sge-slurm-conversion for SGE to SLURM conversion)
            cat ("#!/bin/bash\n",
                 "#SBATCH -o ", scratch.fs, prefix, "-dump/corrcalc-out-%a.txt    # std out\n",
                 "#SBATCH --export=ALL                                   # propagate enrivonment\n",
                 "#SBATCH -n 1                                           # tasks, 1 for R\n",
                 "#SBATCH -p ", cluster.queue, "                         # partition/queue to use\n",
                 "#SBATCH -t 48:0:0                                      # time (2d), required on some systems\n",
                 "#SBATCH --mem=8G                                       # memory \n",
                 "#SBATCH -a 1-", LSF.mut.jid.max, "                     # job array\n",
                 "#SBATCH -J ", prefix, "corr", "                        # job array name\n\n",
                 file=script, sep='')
            cat ("Rscript cna-analysis.r $SLURM_ARRAY_TASK_ID", LSF.mut.jid.max, prefix, mrna, cna, pome, "\n",
                 file=script, append=TRUE)  
          },
          {
            # default -- run serially; may crash or fail
            cat ("#!/bin/bash\n", file=script, sep='')
            cat ("Rscript cna-analysis.r 1 1", prefix, mrna, cna, pome, "\n", file=script, append=TRUE)  
          }
  )
}


run.bsub.file <- function (prefix, force=FALSE) {
  # if force==TRUE, correlation calcs run even if result directories exist
  # if force==FALSE, run only if result directories don't exist
  
  dump.dir <- paste (scratch.fs, prefix, '-dump', sep='')
  out.dir <- paste (scratch.fs, prefix, '-output', sep='')
  run.cmd <- force
  # run if results directories don't exist
  if (any (!file.exists (c (dump.dir, out.dir)))) run.cmd <- TRUE

  if (run.cmd) {
    # remove dump and output directories and create new ones
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
  } else {
    cat ('Using exisiting results for', prefix, '... \n')
  }
}


create.subset <- function (data, sub, fname, na.max=0.5) {
  # create data subset; ensure there aren't too many NAs after subsetting
  data.sub <- data [, c (TRUE, sub)]
  nas <- apply (data.sub[-1], 1, function (x) sum (is.na (x))) / (ncol (data.sub) - 1)
  data.sub <- data.sub [nas < na.max, ]
  write.csv (data.sub, fname, row.names=FALSE)
  return (fname)
}


run.corr.calcs <- function (cls=NULL, groups=NULL, mrna.file=file.path (netgestalt.dir, 'mrna-matrix.csv'),
                            cna.file=file.path (netgestalt.dir, 'cna-matrix.csv'),
                            pome.file=file.path (netgestalt.dir, paste (type, '-matrix.csv', sep='')),
                            force=FALSE) {
  
  mrna.data <- read.csv (mrna.file)
  cna.data <- read.csv (cna.file)
  pome.data <- read.csv (pome.file)
  
  # defaults for cls and groups
  if (is.null (groups)) groups <- 'all'
  if (is.null (cls)) cls <- rep ('all', ncol(pome.data)-1)
  
  for (x in groups) {
    subsamp <- cls == x
    if (sum (subsamp) >= min.cna.N) {
      # run only if there are enough samples
      write.bsub.file (x, 
                       create.subset (mrna.data, subsamp, paste (x, '-mrna-matrix.csv', sep='')),
                       create.subset (cna.data, subsamp, paste (x, '-cna-matrix.csv', sep='')),
                       create.subset (pome.data, subsamp, paste (x, '-pome-matrix.csv', sep=''))
                       )
      run.bsub.file (x, force=force)
    } else {
      cat ('Not enough samples for CNA analysis: Skipping', x, '... \n')
    }
  }
}


finish.plots <- function (groups=NULL) { 
  # wait for all jobs to finish and then run the plotting program
  if (is.null(groups)) groups <- 'all'
  
  jobs.done <-FALSE
  while (! jobs.done) {
    jobs.done <- all ( sapply (paste (groups, 'corr', sep=''), compute.job.done) )
    Sys.sleep (120)
  }
  
  for (x in groups) {
    if (file.exists (script.name (x))) {
      cmd <- switch (compute.cluster.type,
                     uger={paste ('qsub -q long -o uger.out -cwd -j y -l h_vmem=16g -V runR-uger.sh Rscript cna-analysis.r 0 ', 
                                  length (list.files (paste (x, '-output', sep='')))/4, ' ', x, ' NULL NULL NULL', sep='')},
                     slurm={paste ('sbatch -p ', cluster.queue, ' -o slurm.out --mem=16G --export=ALL runR-slurm.sh Rscript cna-analysis.r 0 ', 
                                   length (list.files (paste (x, '-output', sep='')))/4, ' ', x, ' NULL NULL NULL', sep='')},
                     {paste ('Rscript cna-analysis.r 0 1 ', x, ' NULL NULL NULL', sep='')}
      )
      system.x (cmd)
    } else {
      # group not run -- skip plots
      cat ('Skipping CNA plots for', x, '... \n')
    }
  }  
}

