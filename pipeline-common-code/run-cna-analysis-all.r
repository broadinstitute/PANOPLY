
source ('run-cna-analysis.r')


## use this file to set up CNA analysis by subgroups/subsets
## default setup uses all samples


## global -- all samples
# start all correlation calc processess
run.corr.calcs (force=TRUE)

# wait for jobs to finish and run plotting program
finish.plots ()


disease.specific.cna.code <- sprintf ('run-cna-analysis-%s.r', disease)
if (file.exists (disease.specific.cna.code))
  source (disease.specific.cna.code)


