# Estimating prevalence with different levels of sample pooling

# Require package for estimation of true prevalence from estimated

############### FUNCTIONS ######################


# Method for splitting sample x into chunks of size n - Use for pooling
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# Take a single sample from the data and calculate p_bar using the pooled prevalence method
calculate_p_bar <- function(data, num_samples, poolsize, Se, Sp){
  this_replicate <- sample(x = data,size = num_samples)
  return( pooled_prevalence(this_replicate, poolsize, Se, Sp) )
}

# Calculate p_bar for a single replicate. First, divide sample into pools. Then, decide if pools are positive or negative.
# Finally, use formula to calculate p_bar
# Here, Se and Sp are defaults, but can also be set as distributions (e.g. uniform)
pooled_prevalence <- function(sample,poolsize,Se=0.95,Sp=1.0){
  # Estimate p_bar from sample
  # Each positive sample has Se chance of being 1, otherwise 0
  # Each negative sample has Sp chance of being 0, otherwise 1
  # Divide sample into poolsize number of brackets
  pooled_samples <- chunk(sample,length(sample)/poolsize)
  pooled_results <- rep(0,length(pooled_samples))
  # For each pool - Apply sensitivity
  for (p in 1:length(pooled_samples)) {
    # Se chance to be positive if 1 is present
    if (1 %in% unlist(pooled_samples[p])) {pooled_results[p] <- rbinom(1,size=1,prob=Se)}
    else {pooled_results[p] <- rbinom(1,size=1,prob=(1-Sp))}
  }
  p_bar <- 1.0 - ((Se - mean(pooled_results))/(Se+Sp-1))^(1/poolsize)
  return(p_bar)
}

# Method for performing the full experiment
# Pseudocode
# For each population (different prevalences)
#   For each number of samples
#     For each level of pooling
#       For each replicate
#         Determine estimated prevalence


run_simulation_experiment.finite <- function(num_samples,prevalence,pooling_levels,replicates,population,Se=0.95,Sp=1.0){
  # Setup progress bar
  i <- 0
  pb <- txtProgressBar(min = 0, max = length(num_samples)*length(prevalence)*length(pooling_levels), style = 3)
  
  # INITIALIZE RESULTS
  Res_parameters <- data.frame(N=vector(),p=vector(),k=vector())
  Res_replicates <- as.data.frame(matrix(ncol=replicates,nrow=0,dimnames=list(NULL,seq(1:replicates))))
  
  for (prev in 1:length(prevalence)){ #1  (nsamp in 1:length(num_samples))
    for (nsamp in 1:length(num_samples)){ #2 (col in 1:ncol(population))
      for (pool in 1:length(pooling_levels)){ #3
        # Calculate p-bar for all replicates of this parameter combination
        this_sample <- replicate(replicates,calculate_p_bar(population[,prev],num_samples[nsamp],pooling_levels[pool],Se,Sp))
        
        # Add parameter combinations
        Res_parameters[nrow(Res_parameters)+1,] <- c(num_samples[nsamp],prevalence[prev], pooling_levels[pool])
        Res_replicates[nrow(Res_replicates)+1,] <- this_sample
        
        i <- i + 1
        setTxtProgressBar(pb, i)
      }
      
    }
  }
  return( cbind(Res_parameters, Res_replicates))
}


accumulate_fractions <- function(replicates,n,p,k,Se=0.95,Sp=1.0){
  # Run a number of replicates and extract low, median and high quantiles from distribution
  # of positive fraction of wells
  my.set <- replicate(replicates,expr = generate_fractions(n,p,k,Se,Sp))
  return(quantile(my.set,probs = c(0.025,0.5,0.975)))
}

generate_fractions <- function(n,p,k,Se=0.95,Sp=1.0){
  # Generating positive and negative pools according to a binomial function
  # Only for infinite populations
  pools <- rbinom(n = n/k,prob = pool_positive(p,k,Se,Sp),size=1)
  fracs <- mean(pools)
  return(fracs)
}

pool_positive <- function(p,k,Se=0.95,Sp=1.0){
  # Function for generating probabilities of a pool testing positive
  # This can happen both from true positives (at least one positive sample in pool)
  # or false positives, ie there are no true positive samples in the pool but the 
  # test says otherwise
  atleastonepos <- (1-(1-p)^k)*Se
  nopos <- (1-p)^k * (1-Sp)
  return(atleastonepos + nopos)
}

tu.method <- function(frac,k,Se=0.95,Sp=1.0){
  # Tu method of finding p-bar, given a fraction of positive pools
  return( 1 - ((Se-frac)/(Se+Sp-1))^(1/k) )
}

run_simulation_experiment.infinite <- function(num_samples,prevalence,pooling_levels,replicates,Se=0.95,Sp=0.99){
  # Setup progress bar
  i <- 0
  pb <- txtProgressBar(min = 0, max = length(num_samples)*length(prevalence)*length(pooling_levels), style = 3)
  
  # INITIALIZE RESULTS
  Res <- data.frame(N=vector(),p=vector(),k=vector(),low=vector(),median=vector(),high=vector())

  for (prev in 1:length(prevalence)){
    for (nsamp in 1:length(num_samples)){
      for (pool in 1:length(pooling_levels)){
        # Calculate p-bar for all replicates of this parameter combination
        # First, calculate the fractions of positive pools
        fractions <- accumulate_fractions(replicates=replicates,p = prevalence[prev],n = num_samples[nsamp],k=pooling_levels[pool],Se,Sp)
        this_sample <- tu.method(fractions,pooling_levels[pool],Se,Sp)
        
        # Add parameter combinations
        Res[nrow(Res)+1,] <- c(num_samples[nsamp],prevalence[prev], pooling_levels[pool],this_sample)
        
        i <- i + 1
        setTxtProgressBar(pb, i)
      }
    }
  }
  return(Res)
}

#### FREEDOM FROM DISEASE ####

sample_size_machine.k <- function(p, k,alpha=0.05,Sp=1.0,Se=0.95){
  # Freedom from disease. How many samples are needed to have a (1-alpha)
  # probability of sampling at least one positive sample, depending on true prevalence?
  #
  # In other words, what is the highest n where we would reasonably expect all tests to be negative
  # All tests need to be negative, ie, true neg or false neg
  falseneg <- (1-Se)*(1-(1-p)^k)  # 1 - Se*(1-(1-p)^k)
  trueneg <- (1-p)^k * (Sp)
  return(ceiling(log(1-(1-alpha))/log(falseneg+trueneg)))
}


# Difference between two binomial functions - Credit https://gist.github.com/coppeliaMLA/9681819
diffBin<-function(z, n1, p1, n2, p2){
  # calculates probability for z (x-y)
  prob<-0
  if (z>=0){  
    for (i in 1:n1){     
      prob<-prob+dbinom(i+z, n1, p1)*dbinom(i, n2, p2)
    }
  }
  else
  {
    for (i in 1:n2){     
      prob<-prob+dbinom(i+z, n1, p1)*dbinom(i, n2, p2)
    }
  }
  return(prob)
}


#### TESTING EFFICIENCY OF VARIOUS STRATEGIES ####

binary_search_fixedsample <- function(sample,k,Se=0.95,Sp=0.99,correct=0,incorrect=0,num_tests=0){
  # Divide sample into sizes of k. Aggregate binary search numbers
  if (length(sample) == k) { splitssamples <- split(sample,f = 1)}
  else {splitssamples <- chunk(sample,length(sample)/k)}
  res <- list("num"=0,"corr"=0,"inc"=0)
  for (s in splitssamples){
    tmp <- binary_search_sample(s,Se=Se,Sp=Sp)
    res$num <- res$num + tmp$num
    res$corr <- res$corr + tmp$corr
    res$inc <- res$inc + tmp$inc
  }
  return(res)
}

binary_search_sample <- function(sample,Se=0.95,Sp=0.99,correct=0,incorrect=0,num_tests=0){
  # Takes in a sample pool of N people. Tests if sample positive or negative
  # If negative - Clear entire sample (and evaluate num_tests, correct and incorrect diagnoses)
  # If positive - Split sample into 2, add up num_tests,correct and incorrect diagnoses
  # num_tests should know the number of tests in this group WITH lower groups (not higher)
  if (1 %in% sample){
    # True sick in sample pool
    testres <- rbinom(n=1,size=1,prob=Se)
    num_tests <- 1
  }
  else {
    # No sick in sample pool
    testres <- rbinom(n = 1,size = 1,prob = 1-Sp)
    num_tests <- 1
  }
  if (testres == 1) {
    # Positive test result
    if (length(sample) == 1){
      # Single sample
      correct <- sum(sample==1)
      incorrect <- sum(sample==0)
      return(list("num"=num_tests,"corr"=correct,"inc"=incorrect))
    }
    else {
      # Binary split
      sample_div <- chunk(x = sample,n = 2)
      left <- binary_search_sample(sample_div[[1]],Se,Sp)
      right <- binary_search_sample(sample_div[[2]],Se,Sp)
      num_tests <- num_tests + left$num + right$num
      correct <- left$corr + right$corr
      incorrect <- left$inc + right$inc
      return(list("num"=num_tests,"corr"=correct,"inc"=incorrect))
    }
  }
  else {
    # Negative test result
    correct <- sum(sample==0)
    incorrect <- sum(sample==1)
    return(list("num"=num_tests,"corr"=correct,"inc"=incorrect))
  }
}

#### PLOTTING ####

# Takes a bunch of rows and draws polygons
generate_polygon_rows <- function(rows){
  poolsizes <- rows$k
  xcoords <- c(poolsizes,rev(poolsizes))
  ycoords <- rep(0,length(xcoords))
  low_estimates <- rows$low
  high_estimates <- rows$high
  for( i in 1:nrow(rows)){
    ycoords[i] <- low_estimates[i]
    ycoords[length(ycoords)-i+1] <- high_estimates[i]
  }
  return(matrix(cbind(xcoords,ycoords),ncol=2))
  
}

# Takes a bunch of rows and draws the actual prevalence and the 2.5% and 97.5% estimates of p-bar
# Only for results from sampling (finite) population
create_plot_row.finite <- function(rows){
  poolsizes <- unique(rows$k)
  prevalence <- unique(rows$p)
  low_estimates <- apply(X=rows[,4:ncol(rows)],FUN=quantile,MARGIN=1,probs = 0.025,na.rm=TRUE)
  high_estimates <- apply(X=rows[,4:ncol(rows)],FUN=quantile,MARGIN=1,probs = 0.975,na.rm=TRUE)

  plot(poolsizes,rep(prevalence,length(poolsizes)),type="l",col="black",lty=2,ylab="Prevalence", xlab="Pooling (k)",xaxt="n",ylim=c(min(low_estimates),max(high_estimates)))
  axis(side=1,at=poolsizes)
  samplesizes <- unique(rows$N)
  mycols <- rainbow(length(samplesizes))
  legend("topright",legend=samplesizes,fill = mycols)
  for( n in 1:length(samplesizes)){
    polygon(generate_polygon_rows(rows[rows$N == samplesizes[n],]),density=0,col=mycols[n])
  }
}

create_plot_row <- function(rows){
  orig.samplesizes <- unique(rows$N)
  mycols <- rainbow(length(orig.samplesizes))
  
  # REMOVE NaN rows
  rows <- rows[!is.na(rows$high),]
  # REMOVE ROWS WHERE ESTIMATES ARE TOO BIG
  rows <- rows[!rows$high > 0.3,]
  poolsizes <- unique(rows$k)
  prevalence <- unique(rows$p)
  # No estimates should be lower than 0
  rows$low <- sapply(X = rows$low,FUN=max,0)
  rows$median <- sapply(X=rows$median,FUN=max,0)
  rows$high <- sapply(X=rows$high,FUN=max,0)
  low_estimates <- rows$low
  mediab_estimtes <- rows$median
  high_estimates <- rows$high
  
  plot(poolsizes,rep(prevalence,length(poolsizes)),type="l",col="black",lty=2,ylab="Prevalence", xlab="Pooling (k)",xaxt="n",log="x",ylim=c(min(low_estimates),max(high_estimates)))
  axis(side=1,at=poolsizes)
  samplesizes <- unique(rows$N)
  mycol.indicator <- orig.samplesizes %in% samplesizes
  mycols <- mycols[mycol.indicator]
  legend("topright",legend=samplesizes,fill = mycols)
  for( n in 1:length(samplesizes)){
    polygon(generate_polygon_rows(rows[rows$N == samplesizes[n],]),density=0,col=mycols[n])
  }
}

# Create plots for all different prevalences
create_plot_TOT <- function(results){
  # Find unique prevalences
  for( prev in unique(results$p)){
    create_plot_row(results[results$p == prev,])
  }
}


#####################################
