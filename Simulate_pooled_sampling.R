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
    else {pooled_results[p] <- 0} # ONLY WHEN SPECIFICITY = 1.0
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


run_simulation_experiment <- function(num_samples,prevalence,pooling_levels,replicates,population,Se=0.95,Sp=1.0){
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
        
        # Add meta-info to this_sample
        cname <- paste0(num_samples[nsamp],"_",prevalence[col],"_",pool)
        
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



# PLOTTING

# Takes a bunch of rows and draws polygons
generate_polygon_rows <- function(rows){
  poolsizes <- rows$k
  xcoords <- c(poolsizes,rev(poolsizes))
  ycoords <- rep(0,length(xcoords))
  low_estimates <- apply(X=rows[,4:ncol(rows)],FUN=quantile,MARGIN=1,probs = 0.025,na.rm=TRUE)
  high_estimates <- apply(X=rows[,4:ncol(rows)],FUN=quantile,MARGIN=1,probs = 0.975,na.rm=TRUE)
  for( i in 1:nrow(rows)){
    ycoords[i] <- low_estimates[i]
    ycoords[length(ycoords)-i+1] <- high_estimates[i]
  }
  return(matrix(cbind(xcoords,ycoords),ncol=2))
  
}

# Takes a bunch of rows and draws the actual prevalence and the 2.5% and 97.5% estimates of p-bar
create_plot_row <- function(rows){
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

# Create plots for all different prevalences
create_plot_TOT <- function(results){
  # Find unique prevalences
  for( prev in unique(results$p)){
    create_plot_row(results[results$p == prev,])
  }
}


#####################################

####### PARAMETERS ##################

# SET PARAMETERS
num_samples <- c(200,500,1000,1500,2000,3000,5000)
prevalence <- c(0.001,0.003,0.01,0.03,0.1)
pooling_levels <- c(1,3,5,7,10,15,20,25,30)
replicates <- 10000
population_size <- 500000
Se <- 0.95 
Sp <- 1.0

# Generate populations (Maybe generate individuals in the tests instead?)
# Population prevalence can start at 0.001, increment 3x, so 0.003, 0.01, 0.03, 0.1
population <- (matrix(ncol=length(prevalence),nrow=population_size))
for ( p in 1:length(prevalence)) {
  pop <- rbinom(n = population_size,size = 1,prob = prevalence[p])
  population[,p] <- pop
}
population <- as.data.frame(population,dimnames=list(NULL,prevalence))

#####################################



#########  SIMULATION  ###############

# PERFORM EXPERIMENT
RES <- run_simulation_experiment(num_samples,prevalence,pooling_levels,replicates,population)

######################################



########### PLOTS ####################

# Plotting all
create_plot_TOT(RES)

# Plotting one prevalence
create_plot_row(RES[RES$p==0.001,])

# Plotting manually
plot(poolsizes,rep(0.03,length(poolsizes)),type="l",col="black",lty=2,ylab="Prevalence", xlab="Pooling (k)",xaxt="n",ylim=c(0,0.07))
axis(side=1,at=poolsizes)
legend("topright",legend=samplesizes,fill = mycols)
mypop <- RES[RES$p==0.1,]
polygon(generate_polygon_rows(mypop[mypop$N==200 & mypop$k < 13,]),density=0,col=mycols[1])
polygon(generate_polygon_rows(mypop[mypop$N==500 & mypop$k < 20,]),density=0,col=mycols[2])
polygon(generate_polygon_rows(mypop[mypop$N==1000 & mypop$k < 20,]),density=0,col=mycols[3])
polygon(generate_polygon_rows(mypop[mypop$N==1500 & mypop$k < 25,]),density=0,col=mycols[4])
polygon(generate_polygon_rows(mypop[mypop$N==2000 & mypop$k < 25,]),density=0,col=mycols[5])
polygon(generate_polygon_rows(mypop[mypop$N==3000 & mypop$k < 30,]),density=0,col=mycols[6])
polygon(generate_polygon_rows(mypop[mypop$N==5000 & mypop$k < 30,]),density=0,col=mycols[7])

