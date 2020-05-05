

####### PARAMETERS ##################

# SET PARAMETERS
num_samples <- c(200,500,1000,1500,2000,3000,5000)
prevalence <- c(0.001,0.003,0.01,0.03,0.1)
pooling_levels <- c(1,3,5,7,10,15,20,25,30,40,50,70,100,200)
replicates <- 10000
population_size <- 500000
Se <- 0.95 
Sp <- 0.99

# Generate populations
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

RES.TU <- run_simulation_experiment.2(num_samples = num_samples,prevalence = prevalence,pooling_levels = pooling_levels,replicates = 10000)

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



##### FREEDOM FROM DISEASE #######



poolsizes <- c(1,3,5,7,10,13,15,20,25,30)

wideprevalence <- seq(0.0001,0.1,0.0001)

my.df <- sapply(X=poolsizes,FUN = sample_size_machine.k,p=wideprevalence,Sp=0.99)

yticks <- c(1,2,3,4,5,7,10,20,50,100,200,500,1000)

plot(0.1,0.1,xlim=c(0,0.1),ylim=c(1,1000),type="n",yaxt="n",xaxt="n",xlab="Prevalence",ylab="Sample size",log="y")
axis(side=1,at=seq(0,0.1,0.01))
axis(side=2,at=yticks)
mycolswide <- rainbow(10)
lines(wideprevalence,my.df[,1],type="l",col=mycolswide[1])
lines(wideprevalence,my.df[,2],type="l",col=mycolswide[2])
lines(wideprevalence,my.df[,3],type="l",col=mycolswide[3])
lines(wideprevalence,my.df[,4],type="l",col=mycolswide[4])
lines(wideprevalence,my.df[,5],type="l",col=mycolswide[5])
lines(wideprevalence,my.df[,6],type="l",col=mycolswide[6])
lines(wideprevalence,my.df[,7],type="l",col=mycolswide[7])
lines(wideprevalence,my.df[,8],type="l",col=mycolswide[8])
lines(wideprevalence,my.df[,9],type="l",col=mycolswide[9])
lines(wideprevalence,my.df[,10],type="l",col=mycolswide[10])
legend("topright",legend=poolsizes,fill=mycolswide,title = "k")
