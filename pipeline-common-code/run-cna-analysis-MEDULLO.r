
source ('run-cna-analysis.r')


## subtypes
sub.cls <- read.cls ( file.path (netgestalt.dir, 'abSubgroup.cls') )
sub.groups <- c ('G3a', 'G3b', 'G4', 'SHHa', 'SHHb')

gr34.cls <- as.vector (sapply (sub.cls, function (x) ifelse (x %in% c ('G3a','G3b', 'G4'), "G34", "non-G34")))
gr34.groups <- c ('G34')

shhb.g4.cls <- as.vector (sapply (sub.cls, function (x) ifelse (x %in% c ('G4','SSHb'), 'G4SHHb', "non-G4SHHb")))
shhb.g4.groups <- c ('G4SHHb')


## major subtypes
msub.cls <- read.cls ( file.path (netgestalt.dir, 'Subgroup.cls') )
msub.groups <- c ('GR3', 'SHH')



# start all correlation calc processess
run.corr.calcs (sub.cls, sub.groups)
run.corr.calcs (gr34.cls, gr34.groups)
run.corr.calcs (shhb.g4.cls, shhb.g4.groups)
run.corr.calcs (msub.cls, msub.groups)



# wait for jobs to finish and run plotting program
finish.plots (sub.groups)
finish.plots (gr34.groups)
finish.plots (shhb.g4.groups)
finish.plots (msub.groups)


