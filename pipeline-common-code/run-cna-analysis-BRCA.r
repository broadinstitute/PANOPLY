
source ('run-cna-analysis.r')


## PAM-50 subtypes
pam50.cls <- read.cls ( file.path (norm.dir, paste (type, subset.str, '-pam50.cls', sep='')) )
keep <- pam50.cls == "Basal" | pam50.cls == "Her2" | pam50.cls == "LumA" | pam50.cls == "LumB"
pam50.cls <- pam50.cls [keep]
pam50.groups <- c ('Basal', 'Her2', 'LumA', 'LumB')



## luminal A+B
lum.cls <- as.vector (sapply (pam50.cls, function (x) ifelse (grepl ('Lum', x), "Luminal", "non-Luminal")))
lum.groups <- c ('Luminal')


## ER+ and ER-
er.cls <- read.cls ( file.path (norm.dir, paste (type, subset.str, '-er.cls', sep='')) )
er.cls <- sapply (er.cls, tolower)
er <- er.cls [ keep ]   # remove non-PAM50 samples
er <- sapply (er, function (x) ifelse (x == 'positive', 'ER-pos', 'ER-neg'))
er.groups <- c ('ER-pos', 'ER-neg')


## PR+ and PR-
pr.cls <- read.cls ( file.path (norm.dir, paste (type, subset.str, '-pr.cls', sep='')) )
pr.cls <- sapply (pr.cls, tolower)
pr <- pr.cls [ keep ]   # remove non-PAM50 samples
pr <- sapply (pr, function (x) ifelse (x == 'positive', 'PR-pos', 'PR-neg'))
pr.groups <- c ('PR-pos', 'PR-neg')


## HER2+ and Triple-Negative
her2.cls <- read.cls ( file.path (norm.dir, paste (type, subset.str, '-her2.cls', sep='')) )
her2.cls <- sapply (her2.cls, tolower)
her2 <- her2.cls [ keep ]   # remove non-PAM50 samples
her2 <- sapply (her2, function (x) ifelse (x == 'positive', 'Her2-pos', 'Her2-neg'))

tn.cls <- er.cls=="negative" & pr.cls=="negative" & her2.cls=="negative"
tn <- tn.cls [ keep ]
for (i in 1:length(tn)) if (tn[i]) her2[i] <- "Triple-neg"

her2.groups <- c ('Her2-pos', 'Triple-neg')


# start all correlation calc processess
run.corr.calcs (pam50.cls, pam50.groups)
run.corr.calcs (er, er.groups)
run.corr.calcs (pr, pr.groups)
run.corr.calcs (her2, her2.groups)
run.corr.calcs (lum.cls, lum.groups)



# wait for jobs to finish and run plotting program
finish.plots (pam50.groups)
finish.plots (er.groups)
finish.plots (pr.groups)
finish.plots (her2.groups)
finish.plots (lum.groups)


