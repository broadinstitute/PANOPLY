###################################################################################################################
## Filename: modT.r
## Created: June 16, 2016
## Authors: Karsten Krug, DR Mani
##
## Purpose: Make functions written by D. R. Mani accessible via Shiny.
##
## The file contains function to perform moderated T and F-test statistics, 2-component normlization
## ratio data and reproducibility filtering (Blandt-Altlman and and a generalized version). All functions
## were implemented by DR Mani and modified by Karsten Krug to be compatible with the Shiny framework.
##
###################################################################################################################


#####################################################################################
##
##                        two sample moderated t-test
##
## - code written by dr mani
## - code modified by karsten krug
##   20151211 'label'
#####################################################################################
modT.test.2class <- function (d, output.prefix, groups, id.col=NULL, data.col=NULL,
                              group.na.rm=FALSE, nastrings=c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?"), label=NULL) {

    cat('\n-- modT.test.2class --\n')

    ## store group names
    groups.org <- groups
    ##groups <- as.numeric(as.factor(groups))
    groups <- as.numeric(factor(groups, levels=sort(unique(groups)))) ## kk 20170106

    ##cat(groups.org, '\n', groups, '\n')

    id <- d[ , id.col]

    ## extract data columns
    if (is.null (data.col)) data <- d [, setdiff (colnames (d), id.col)]
    else data <- d [, make.names (data.col)]


    ## moderated t test for 2 classes
    design.mat <- cbind (ref=1, comparison=groups)
    mod.t.result <- moderated.t (data, design.mat)

    ## 20151211 kk
    mod.t.result <- data.frame( mod.t.result, Log.P.Value=-10*log(mod.t.result$P.Value,10), stringsAsFactors=F)

    ## add label
    if(!is.null(label))
        colnames(mod.t.result) <- paste(colnames(mod.t.result), label, sep='.')

    mod.t <- data.frame ( cbind (data.frame (id=id), data, mod.t.result), stringsAsFactors=F )
    ##mod.t <- data.frame ( cbind (data.frame (id),  mod.t.result, data) )
    rownames(mod.t) <- id ##make.unique(as.character(mod.t[, 1]), sep='_')
    ##colnames(mod.t)[1] <- 'id'

    ##write.csv (final.results, paste (output.prefix, ".csv", sep=''), row.names=FALSE)

    ## write out / return results
    final.results <- mod.t

    ##invisible (final.results)
    return( list(input=d, output=final.results, groups=groups.org) )
}


######################################################################################################
##                               One-sample moderated t-test
##
##
## run moderated t-test, and plot results
## mainly for iTRAQ, but can be used of other data
##
## code written by mani dr
## code modified by karsten krug
##
## 20151210 parameter 'label'
######################################################################################################
modT.test <- function (d, output.prefix, id.col=NULL, data.col=NULL, fix.id=FALSE,
                       p.value.alpha=0.05, use.adj.pvalue=TRUE, apply.log=FALSE,
                       na.rm=FALSE, nastrings=c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?"),
                       plot=TRUE, pairs.plot.2rep=FALSE, limits=NULL, xlab="", ylab="", label='', ...) {
  #
  # data.file should contain one peptide in each row.
  # The columns contain the normalized log-ratio from each replicate
  # (technical or biological). The ratio is based on classes of interest
  # that need to be distinguished: i.e., ratio = intensity_A / intensity_B;
  # this test calculates the p-value for determining if peptide p is
  # differentially regulated between classes A and B (i.e., if the log
  # ratio is different from 0).
  # While the standard scatter plot routine can only handle 2 replicates,
  # a pairs plot is created when there are more than 2 replicates.
  # The moderated t-test can be applied to any number of replicates.
  # An id column can be optionally included in the data.file to track
  # peptides (row numbers are used as id if a column is not specified).
  #
  # graphics can be controlled using ...
  #  when using scatterhist (for 2 replicates), this can be arguments to points
  #  when > 2 replicates are present, ... can include arguments to points in
  #   addition to: plot.col, subset.col, hist.col, hist.breaks,
  #                prefix (for correlation), cex.cor

     cat('\n-- modT.test --\n')
  id <- d[, id.col]
  ##if ( any (duplicated (id)) ) stop ('IDs are not unique. Use fix.id=TRUE option')

  # extract data columns
  if (is.null (data.col)) data <- d [, setdiff (colnames (d), id.col)]
  else data <- d [, make.names (data.col)]

    data <- data.matrix(data)

  # log transform is required
  if (apply.log) data <- log2 (data)

    ## moderated t test
      ##View(data)
    mod.t.result <- moderated.t (data)
    ##View(data)
  if (use.adj.pvalue) mod.sig <- mod.t.result [,'adj.P.Val'] <= p.value.alpha
  else  mod.sig <- mod.t.result [,'P.Value'] <= p.value.alpha
  change <- apply (data, 1,
                   function (x) {
                     x <- x [is.finite (x)]
                     ret.value <- '?'
                     if ( all (x < 0) ) ret.value <- 'down'
                     else if ( all (x > 0)) ret.value <- 'up'
                     return (ret.value)
                   })
    ## 20151210 kk
    mod.t.result <- data.frame( mod.t.result, change=change, significant=mod.sig, Log.P.Value=-10*log(mod.t.result$P.Value,10), stringsAsFactors=F)

    ## add label
    if(!is.null(label))
        colnames(mod.t.result) <- paste(colnames(mod.t.result), label, sep='.')

    ##mod.t <- data.frame ( cbind (data.frame (id), data, mod.t.result, change=change, significant=mod.sig) )
    mod.t <- data.frame ( cbind (data.frame (id=id), data, mod.t.result), stringsAsFactors=F )
    ##colnames (mod.t)[1] <- id.col   # retain id.col (if provided)
    ##rownames(mod.t) <- make.unique( as.character(mod.t[,1]), sep='_' )
    rownames(mod.t) <- id
   ## colnames(mod.t)[1] <- 'id'

    final.results <- mod.t
    cat('\n-- modT.test exit --/n')
    return( list(input=d, output=final.results) )
}


#################################################################################################################
## run moderated F-test, and plot results
## mainly for iTRAQ, but can be used for other data
##
## code written by mani dr
## code modified by karsten krug
modF.test <- function (d, class.vector, output.prefix, id.col=NULL,
                       p.value.alpha=0.05, use.adj.pvalue=TRUE,
                       na.rm=FALSE, nastrings=c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?"),
                       plot=TRUE, limits=NULL, xlab="", ylab="", plot.by.group=TRUE,
                       add.xy.axes=TRUE, ...) {
  #
  # data.file should contain one peptide in each row.
  # The columns contain the normalized log-ratio from each replicate
  # (technical or biological), representing a group/class/comparison of
  # interest. The ratio is based one class of interest and its reference
  # or control i.e., ratio = intensity_class / intensity_ref;
  # this test calculates the p-value for determining if peptide p is
  # differentially regulated in any of the classes (i.e., if the log
  # ratio is different from 0 for any of the classes).
  # The class.vector provides the class label for each non-id column.
  # A pairs plot is generated for sample from each class.
  # An id column can be optionally included in the data.file to track
  # peptides (row numbers are used as id if a column is not specified).
  #
  # graphics can be controlled using ...
  #  when using scatterhist (for 2 replicates), this can be arguments to points
  #  when > 2 replicates are present, ... can include arguments to points in
  #   addition to: plot.col, subset.col, hist.col, hist.breaks,
  #                prefix (for correlation), cex.cor

    id <- id <- d[,id.col]
    data <-  d [, setdiff (colnames (d), id.col)]

    cat('\n-- modF.test --\n')

  # create design matrix
  f <- factor (class.vector)
  design <- model.matrix ( ~ 0 + f )

  # moderated F test
  fit <- lmFit (data, design)
  fit <- eBayes (fit)

  sig <- topTable (fit, number=nrow(data), sort.by='none')
  if (use.adj.pvalue) mod.sig <- sig [,'adj.P.Val'] <= p.value.alpha
  else  mod.sig <- sig [,'P.Value'] <= p.value.alpha
  non.na.n <- apply (data, 1, function (x) { sum (is.finite (x)) })


    mod.f <- data.frame ( cbind (id=id, sig, significant=mod.sig, total.n=non.na.n, Log.P.Value=-10*log(sig[,'P.Value'] ,10)), stringsAsFactors=F )

    ##if(!is.null(label))
    ##    colnames(mod.f) <- paste(colnames(mod.f), label, sep='.')

    final.results <- mod.f
    colnames(final.results) <- sub("^f", "logFC.", colnames(final.results))

    cat('\n-- modF.test exit --\n')
    return( list( input=d, output=final.results, groups=class.vector)  )
  ##invisible (final.results)
}

###########################################################################
##
##           Moderated t-test for significance testing
##
## code written by Mani DR
## code modified by karsten krug
##    - single function for one/two sample test
##
##########################################################################
moderated.t <- function (data, design=NULL) {
    ## data is a table with rows representing peptides/proteins/genes
    ## and columns representing replicates

     cat('\n-- moderated.t --\n')

    data.matrix <- data.frame (data, stringsAsFactors=F)
    ## the design matrix is expected to be:
    ##    ref    comparison
    ##     1         0
    ##    1         0
    ##         ...
    ##     1         1
    ##  where comparison has 0's and 1's corresponding
    ##  to the column position of the 2 classes in data
    ##  (see limma user manual section 13)

    #cat('here1 ')
    ##save(data.matrix, file='data.matrix.RData')
		
		## use robust eBayes
		## see Phipson, B., Lee, S., Majewski, I. J., Alexander, W. S., & Smyth, G. (2013). Tech Report.
		##     Empirical Bayes in the presence of exceptional cases, with application to microarray data.
    
		#############################################
    ## two sample test
    if(!is.null(design)){
        m <- lmFit (data.matrix, design)
        m <- eBayes (m, robust=TRUE)
        sig <- topTable (m, coef=colnames (design)[2], number=nrow(data), sort.by='none')
    } else {
    #############################################
    ## one sample test
        ##cat('here2 ')
        m <- lmFit (data.matrix, method='robust')
        ##cat('here3 ')
        m <- eBayes (m, robust=TRUE)
        ##at('here4 ')
        sig <- topTable (m, number=nrow(data), sort.by='none')
        ##cat('here5 ')
    }
     cat('\n-- moderated.t exit --\n')
  return (sig)
}

############################################################################################
##
##              Generalized reprodicibility filter for > 2 replicates
##
## written by D R Mani
############################################################################################
reproducibility.filter <- function (data, id.col='id', alpha=0.05) {
  ##
  ## Reproducibility Filter using lme4
  ## Theory: MethComp book (pp 58-61). Comparing Clinical Measurement Methods by Bendix Carstensen
  ## Implementation: MethComp book pg 142, but don't include item (=id) in the fixed effects
  ## -- this is unnecessary for the application and makes computations very time consuming;
  ## all we really need to assess reproducibility are the variances
  ##
  ## NB: using library (nlme) and lmer is much more convoluted since incorporating the
  #      var-cov matrix stratified by replicate (=method) is not easy
  #      (see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2007q3/000248.html)
  #      something like:
  #        model <- lmer (y ~ rep + (rep1|id)+(rep2|id)+..., data=data.long)
  #      is possible, but the above does not work for 2 replicates (gives same results for >2 reps)
  #
    cat('\n-- reproducibility.filter --\n')

  d <- data [, setdiff (colnames (data), id.col)]     # data part of input
  data.long <- melt (data.frame (data, stringsAsFactors=F), id=id.col)    # convert to long format
  colnames (data.long) <- c ('id', 'rep', 'y')
  # keep column order in data so that (i,j) below correctly corresponds to columns
  data.long [,'rep'] <- factor (data.long[,'rep'], levels=colnames (d))

  # exclude missing data points (only missing measurements are removed instead of entire rows)
  data.long <- data.long [ !is.na (data.long[,'y']), ]

  # Model: y_mi = a_m + c_mi + e_mi,  c_mi ~ N(0,tau_m^2), e_mi ~ N(0, sigma_m^2)
  # where m = method and i = item (=id)
  # [Eq 5.2, pg 58, MethComp book]. Also see interpretation of effect on pg 59-61
  model <- lme (y ~ rep,
                random=list (id=pdIdent(~rep)),
                weights=varIdent(form=~1|rep),
                data=data.long)
  n <- nlevels (data.long[,'rep'])
  p <- length (unique (data.long[,'id']))
  df <- p - 1    # approx df for confidence interval (p=# of independent items)

  rep.all <- rep (TRUE, nrow (d))  # vector summarizing reproducibility of each input id
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # variance of method_i - method_j: pg 58
      # var (y_i0 - y_j0) = tau_i^2 + tau_i^2 + sigma_i^2 + sigma_j^2
      # where tau is the sd of the between-item variation for a method
      # and sigma is the sd of the within-item variation for the method
      tau <- as.numeric (unlist (VarCorr(model)[c(i,j),'StdDev']))  # returns tau_i and tau_j
      # VarCorr(model)[n+1,'StdDev'] is sigma_1 (ie sigma for method 1)
      # for methods 2-n, coef (model$modelStruct$varStruct, uncons=F, allcoef=T) has the
      #  multiplying factor to obtain sigma_m using sigma_1
      sigma <- as.numeric (unlist (VarCorr(model)[n+1,'StdDev']))   #
      sigma1 <- sigma * ifelse (i==1, 1, coef (model$modelStruct$varStruct, uncons=F, allcoef=T)[i-1])
      sigma2 <- sigma * coef (model$modelStruct$varStruct, uncons=F, allcoef=T)[j-1]
      total.sd <- sqrt (tau[1]^2 + tau[2]^2 + sigma1^2 + sigma2^2)

      # bias of method_i - method_j: alpha_i - alpha_j
      alpha1 <- ifelse (i==1, 0, fixef(model)[i])
      alpha2 <- fixef (model)[j]
      bias <- alpha1 - alpha2

      # limits of agreement (assuming approx df = p-1)
      t.crit <- qt ( (1-alpha/2), df ) * sqrt ( (p+1) / p )
      ci.low <- bias - t.crit * total.sd
      ci.high <- bias + t.crit * total.sd

      # record reproducibility for method_i - method_j
      rep.ij <- (d[,i] - d[,j]) >= ci.low & (d[,i] - d[,j]) <= ci.high
      # if data is missing, assume that data is reproducible
      rep.ij [ is.na (rep.ij) ] <- TRUE
      rep.all <- rep.all & rep.ij

      # print bias and LoA for sanity check
      ##cat ('rep', i, ' - rep', j, ': bias=', bias, ' ci=(', ci.low, ',', ci.high, ')\n', sep='')
    }
  }
    ## karsten krug 20160301
    ## return rownames of data matrix
    cat('\n-- reproducibility.filter exit --\n')
    return(rep.all)
    ##return( rownames(data)[ which(!rep.all) ] )
}

###########################################################################################
##
##                    two-component mixture model normalization
## written by D R Mani
##
##########################################################################################

two.comp.normalize <- function (sample, type='default', mode.lower.bound=-3) {
  #   1. For all sample types, fit a 2-component gaussian mixture model using normalmixEM.
  #   2. For the bimodal samples, find the major mode M1 by kernel density estimation
  #     2a. Fit the model with one component mean constrained to equal M1
  #     2b. Normalize (standardize) samples using mean (M1) and resulting std. dev.
  #   3. For unimodal (default) samples, find the mode M using kernel density estimation
  #     3a. Fit the model with mean for both components constrained to be equal to M
  #     3b. Normalize (standardize) samples using mean M and smaller std. dev. from model fit
  #
  #  the major mode should be located at a value larger than mode.lower.bound 
  
  # WARNING:
  # This code has a lot of hacks to fix the flakiness of normalmixEM, and the idiosyncracies
  # of the actual data. Carefully re-examine code for new or altered input data
  # Currently works for log-ratio data approximately centered around 0
  cat('\n-- two.comp.normalize --\n')
  is.error <- function(x) inherits(x, "try-error")             # function to check for error
  
  data <- sample [ !is.na (sample) ]
  data.range <- diff (range (data))
  dens <- try (density (data, kernel='gaussian', bw='SJ'))     # gaussian kernel with S-J bandwidth
  if (is.error (dens))                                         # sometimes, SJ bw estimation fails
    dens <- density (data, kernel='gaussian', bw='ucv')        # in such cases, use unbiased CV
  # (see Venalbles & Ripley, 2002, pg, 129
  #  and Density Estimation, S.J.Sheather, Stat. Sci. 2004)
  # find major (highest) mode > -3 (to avoid problems with lower mode having higher density than higher mode)
  x.range <- dens$x > mode.lower.bound  
  dens.x <- dens$x [x.range];  dens.y <- dens$y [x.range]
  mode <- dens.x[which.max(dens.y)]                                                        
  if (type=='bimodal') mean.constr <- c (NA, mode) else mean.constr <- c (mode, mode)
  model <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
  model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
  model.alt <- Mclust (data, G=2, modelNames=c ("V","E"))
  # V results is separate SDs for each cluster; E fits a single SD for both clusters
  if (length (model.alt$parameters$variance$sigmasq)==1)  # latter code expects two SD values
    model.alt$parameters$variance$sigmasq <- rep (model.alt$parameters$variance$sigmasq, 2)
  alt.mu <- model.alt$parameters$mean
  alt.sd <- sqrt (model.alt$parameters$variance$sigmasq)
  # find reproducible model fit that is close to Mclust fit
  # if not, re-fit model -- without this condition
  # normalmixEM produces one-off model fits 
  n.try <- 1
  if (type=='bimodal') model.mode <- which(model$mu==mode)
  else model.mode <- which(model$mu==mode)[which.min (model$sigma)]
  model.other <- model.mode %% 2 + 1
  alt.mode <- ifelse (diff (alt.mu) < data.range*0.05,          # means are close -- 
                      which.min (alt.sd),                       # use min sd to pick alt.mode
                      which.min(abs(model.alt$par$mean-mode)))  #  else use alt.mu closest to mode
  # always using latter can result in consistently picking the wrong alt.mode => no convergence
  alt.other <- alt.mode %% 2 + 1
  while ( abs (model$mu[model.mode] - alt.mu[alt.mode]) > data.range*0.05 || 
          abs (model$sigma[model.mode]-alt.sd[alt.mode]) > data.range*0.05 ||
          model$sigma[model.mode] < 0.1 ||
          (type=='bimodal' && (abs (model$mu[model.other] - alt.mu[alt.other]) > data.range*0.25)) ||
          abs (sum (c (model$mu, model$sigma) - c (model.rep$mu, model.rep$sigma))) > 1e-3 ) {
    # if major mode (and SD of mode) is not within 5% of data range, or if the other mean (for bimodals only) 
    # is not within 25% of the Mclust result, try again
    model <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
    model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
    
    if (n.try > 1000) stop (paste ("Can't fit mixture model ... giving up\n"))
    n.try <- n.try + 1
  }
  
  
  if (type=='bimodal') {
    # sometimes (esp. in phosphoproteome) the minor (lower) mode can be larger than the major (higher) mode
    # this situation is not possible in the unimodal samples
    corrected.mode <- model$mu [which.max(model$mu)]
    if (corrected.mode != mode) {
      cat ('  Lower mode larger than higher mode\n')
      mode <- corrected.mode
    }
  }
  norm.mean <- mode
  norm.sd <- ifelse (type=='bimodal', model$sigma[which(model$mu==mode)], min (model$sigma))
  
  # normalize by standardizing
  data <- data - norm.mean
  data <- data / norm.sd
  
  # return normalized data reorganized to original order
  sample [ !is.na (sample) ] <- data
  cat('\n-- two.comp.normalize exit --\n')
  return ( list (norm.sample=sample, norm.mean=norm.mean, norm.sd=norm.sd, fit=unlist (c(model$mu, model$sigma))) )
}

## 20180417
two.comp.normalize.old <- function (sample, type) {
  #   1. For all sample types, fit a 2-component gaussian mixture model using normalmixEM.
  #   2. For the bimodal samples, find the major mode M1 by kernel density estimation
  #     2a. Fit the model with one component mean constrained to equal M1
  #     2b. Normalize (standardize) samples using mean (M1) and resulting std. dev.
  #   3. For unimodal samples, find the mode M using kernel density estimation
  #     3a. Fit the model with mean for both components constrained to be equal to M
  #     3b. Normalize (standardize) samples using mean M and smaller std. dev. from model fit

  # WARNING:
  # This code has a lot of hacks to fix the flakiness of normalmixEM, and the idiosyncracies
  # of the actual data. Carefully re-examine code for new or altered input data
    cat('\n-- two.comp.normalize --\n')

  data <- sample [ !is.na (sample) ]
  dens <- density (data, kernel='gaussian', bw='SJ')     # gaussian kernel with S-J bandwidth
                                                         # (see Venalbles & Ripley, 2002, pg, 129)
  # find major (highest) mode > -3 (to avoid problems with lower mode having higher density than higher mode)
  x.range <- dens$x > -3
  dens.x <- dens$x [x.range];  dens.y <- dens$y [x.range]
  mode <- dens.x[which.max(dens.y)]
  if (type=='bimodal') mean.constr <- c (NA, mode) else mean.constr <- c (mode, mode)
  model <- normalmixEM (data, k=2, mean.constr=mean.constr)
  model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr)
  model.alt <- Mclust (data, G=2, modelNames="V")
  alt.mu <- model.alt$parameters$mean
  alt.sd <- sqrt (model.alt$parameters$variance$sigmasq)
  # find reproducible model fit that is close to Mclust fit
  # if not, re-fit model -- without this condition
  # normalmixEM produces one-off model fits
  n.try <- 1
  if (type=='unimodal') model.mode <- which(model$mu==mode)[which.min (model$sigma)]
  else model.mode <- which(model$mu==mode)
  model.other <- model.mode %% 2 + 1
  alt.mode <- which.min(abs(model.alt$par$mean-mode)); alt.other <- alt.mode %% 2 + 1
  while ( abs (model$mu[model.mode] - alt.mu[alt.mode]) > 3e-1 || abs (model$sigma[model.mode]-alt.sd[alt.mode]) > 3e-1 ||
          model$sigma[model.mode] < 0.1 ||
          (type=='bimodal' && (abs (model$mu[model.other] - alt.mu[alt.other]) > 1e1)) ||
          !all (c (model$mu, model$sigma) - c (model.rep$mu, model.rep$sigma) < 1e-3) ) {
    # if major mode (and SD of mode) is not within 0.3, or if the other mean (for bimodals only)
    # is not within 1 of the Mclust result, try again
    model <- normalmixEM (data, k=2, mean.constr=mean.constr)
    model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr)

      if (n.try > 100){
          return("No_success")
          ##stop (paste ("Can't fit mixture model ... giving up\n"))
      }
    n.try <- n.try + 1
  }


  if (type=='bimodal') {
    # sometimes (esp. in phosphoproteome) the minor (lower) mode can be larger than the major (higher) mode
    # this situation is not possible in the unimodal samples
    corrected.mode <- model$mu [which.max(model$mu)]
    if (corrected.mode != mode) {
      cat ('  Lower mode larger than higher mode\n')
      mode <- corrected.mode
    }
  }
  norm.mean <- mode
  norm.sd <- ifelse (type=='bimodal', model$sigma[which(model$mu==mode)], min (model$sigma))

  # normalize by standardizing
  data <- data - norm.mean
  data <- data / norm.sd

  # return normalized data reorganized to original order
    sample [ !is.na (sample) ] <- data

    cat('\n-- two.comp.normalize exit --\n')

  return ( list (norm.sample=sample, norm.mean=norm.mean, norm.sd=norm.sd, fit=unlist (c(model$mu, model$sigma))) )
}
