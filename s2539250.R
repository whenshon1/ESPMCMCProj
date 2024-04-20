#Will Henshon
#Student number: s2539250

#This code generates two Bayesian models to simulate the effect of the fungicide
#vinclozolin on human androgen receptors, measured in Chinese Hamster ovary
#cells, from 5 experimental runs of a series of 9 dose sizes. These models are
#smooth and monotonic, since we expect the activity of the receptors to decrease
#monotonically with the concentration. To do this, I create latent variables
# (x) for each of the concentration levels with Normal priors, and then trans-
#form them deterministically into mu, which is a smoothly decreasing sequence
#of positive numbers. I use a uniform prior for M, the expected value of
#the effect at 0 dose for each experimental run. I then model the effect as
#Normal, with a mean of M_j * mu_i, where M_j is M at the jth experimental run
#and mu_i is the mu value for the ith dose index. (I also have uniform priors 
#for tau and tau_0, which reflect the precisions in the Normal priors for the x 
#and effect parameters, respectively). The second model differs in that the mu
#parameters are found for each experimental run separately, so they are indexed
#mu_ij; the rest of the model remains the same. To implement the models, I use
#JAGS for Gibbs sampling. I create trace plots for the highest and lowest 
#effective sample size variables for the first model, and then plots for both
#models that include the experimental data as dots, and the model mean and 
#95% credible intervals for each experimental run as lines and shaded regions.
#Finally, I use DIC to compare the models, concluding that the second model is
#more effective.

#Import jags and set the working directory

#pdf("mono.pdf",height=10,width=6)
#par(mfrow=c(2,1))

library(rjags)
library(ggplot2)
setwd('/Users/willhenshon/Desktop/ExtendedStatsSoloWork/Practical5/')

#Read in the data
vin <- read.table('vin.txt')

#Create an index list in terms of experiment numbers 1-5 instead of the actual
#experiment numbers; this is useful for indexing by experiment
experiments <- unique(vin$exper)
experno <- c(1:length(experiments))
expct <- rep(0,length(vin$exper))
for(i in 1:length(experno)){
  ii <- which(vin$exper==experiments[i])
  expct[ii] <- experno[i]
}

#Generate a vector of the concentration values
conc_length <- length(vin$conc)/length(experiments)
conc_data <- vin$con[1:conc_length]

#Reorganize the experimental effect data into a 5x9 matrix so it can be 
#indexed by experiment number and concentration
effect_mat <- matrix(vin$effect,nrow=length(conc_data), ncol=(5))

#Define I and J, representing loop iterations for the concentration index
#and experiment number index, respectively, for use in the JAGS file
I <- conc_length
J <- length(experno)

#Run JAGS to build the model
mod<-jags.model("M1s2539250.jags", data = list(eff=effect_mat, I=I, J=J))

#Create a coda object with the important parameters from the JAGS file. The 
#sample size of 30000 with a burn-in of 7500 avoids issues with mixing in the 
#first 7000 iterations for some of the smaller-sample size parameters.
sam.coda <- coda.samples(mod,c("M","tau0","mu[2:9]","tau",'fx'), n.iter=30000, 
                         n.burnin=7500)

#Find the effective sample sizes of each parameter
eff_size <- effectiveSize(sam.coda)

#Find the parameters with the largest and smallest sample sizes; for the run 
#selected, these were found to be tau0, the precision for the output prior, 
#with an effective sample size of 15147.15, and fx[9,3], the effect output for 
#the 3rd experiment at the 9th dose index, with a sample size of 464.95. This 
#changes from run-to-run because the algorithm involves randomness. I included 
#effect outputs in this analysis because even though they are not explicitly a
#model parameter, since these are the values we are simulating, it is still 
#important to check whether they all meet the requirements for effective sample 
#size and mixing. Of the parameters, the one with the smallest effective sample 
#size was mu[9] (the mu value for the 9th dose index). 
#
max_size <- which(effectiveSize(sam.coda)==max(eff_size))
min_size <- which(effectiveSize(sam.coda)==min(eff_size))

#Generate trace plots for the largest and smallest sample sized elements. Since
#there appears to be fast mixing for both of these, we conclude that the sample
#size and (lack of) burn-in is sufficient for the model.
traceplot(sam.coda[1][,max_size[[1]]],
          main='Maximum Effective Sample Size Trace Plot')
traceplot(sam.coda[1][,min_size[[1]]],
          main='Minimum Effective Sample Size Trace Plot')


#Create a column in the dataset for the concentration index
dose_ind = c(1:45); dose_ind = dose_ind %% 9
dose_ind[which(dose_ind==0)]=9
vin$dose_ind = dose_ind

#Create a column in the dataset to store the experiment number as a categorical
#variable.
experiments = vin$exper
for(i in 1:length(vin$exper)){
  experiments[i]=as.character(vin$exper[i])
}
  
#Find the credible intervals for the sampled parameters
cred_ints <- apply(sam.coda[[1]],2,quantile,prob=(c(0.025,.975)))

#Create matrices to store the expected effect, and the lower and upper bounds
#of the 95% credible interval, for each experiment and dose index.
eff_mod_output <- matrix(0,9,5)
mod_lower_out <- matrix(0,9,5)
mod_upper_out <- matrix(0,9,5)

#Iterate through the dose indices and experiments
for(i in 1:9){
  for(j in 1:5){
      #Add the credible intervals and effect means for each experiment and dose
      #index to the proper matrix location. Since the effects are Normally
      #distributed, we can find the expected value with the center of the
      #lower and upper value of the credible intervals.
      mod_lower_out[i,j] = cred_ints[[11+18*(j-1)+2*(i-1)]]
      mod_upper_out[i,j] = cred_ints[[12+18*(j-1)+2*(i-1)]]
      eff_mod_output[i,j] = (mod_lower_out[i,j]+mod_upper_out[i,j])/2
  }
}


#Transform the matrices into vectors, and store them to the dataset
mod_eff <- c(eff_mod_output)
lower <- c(mod_lower_out)
upper <- c(mod_upper_out)
vin$mod_eff <- mod_eff
vin$lower <- lower
vin$upper <- upper


#Create a plot
a<- ggplot(vin, aes(x=dose_ind, y=effect))
#Add the experimental data as points
a <- a+ geom_point(aes(colour=experiments))
#Add the simulated data as lines
a <- a+ geom_line(aes(x=dose_ind, y = mod_eff, colour = experiments) )
#Add the 95% credible intervals as shaded regions with dashed borders
a <- a +geom_ribbon(aes(ymin = lower, ymax=upper, colour = experiments),
                    linetype=2, alpha=0.2)
#Add labels to the graph and output it
a <- a+labs(x='Dose Index', y='Effect', 
            title = "Data and Model Output for Model #1")
a <- a + scale_x_continuous(name="Dose Index", breaks=c(1,2,3,4,5,6,7,8,9))
print(a)


#Create the second version of the model, and generate a coda file with the
#necessary parameters. I use a burn-in of 20000 for this model because 
#the mixing improves for the smaller effective sample size variables after that
#point.
mod2<-jags.model("M2s2539250.jags", data = list(eff=effect_mat, I=I, J=J))
sam2.coda <- coda.samples(mod2,c("M","tau0","mu[2:9,1:5]","fx"),n.iter=50000, 
                          n.burnin=20000)

#Generate the credible intervals for the parameters
cred_ints2 <- apply(sam2.coda[[1]],2,quantile,prob=(c(0.025,.975)))

#Create matrices to store the model outputs and the credible intervals
eff_mod_output2 <- matrix(0,9,5)
mod_lower_out2 <- matrix(0,9,5)
mod_upper_out2 <- matrix(0,9,5)

#Iterate through the dose indices and experiments
for(i in 1:9){
  for(j in 1:5){
    #Add the credible intervals and effect means for each experiment and dose
    #index to the proper matrix location
    mod_lower_out2[i,j] =cred_ints2[[11+18*(j-1)+2*(i-1)]]
    mod_upper_out2[i,j] =cred_ints2[[12+18*(j-1)+2*(i-1)]]
    eff_mod_output2[i,j] = (mod_lower_out2[i,j] +mod_upper_out2[i,j])/2
  }
}

#Transform the matrices into vectors and store them to the dataset
mod_eff2 <- c(eff_mod_output2)
lower2 <- c(mod_lower_out2)
upper2 <- c(mod_upper_out2)
vin$mod_eff2 <- mod_eff2
vin$lower2 <- lower2
vin$upper2 <- upper2

#Create a plot for the second simulation
a2<- ggplot(vin, aes(x=dose_ind, y=effect))
#Add the experimental data as points
a2 <- a2+ geom_point(aes(colour=experiments))
#Add the simulated data as lines
a2 <- a2+ geom_line(aes(x=dose_ind, y = mod_eff2, colour = experiments) )
#Add the credible intervals as shaded regions
a2 <- a2 +geom_ribbon(aes(ymin = lower2, ymax=upper2, colour = experiments),
                    linetype=2, alpha=0.2)
#Add labels and output the graph
a2 <- a2+labs(y='Effect', 
            title = "Data and Model Output for Model #2")
a2 <- a2 + scale_x_continuous(name="Dose Index", breaks=c(1,2,3,4,5,6,7,8,9))

print(a2)


#Regenerate the models with 2 chains, which is required for the deviances
#to be calculated with DIC
mod_1<-jags.model("M1s2539250.jags", data = list(eff=effect_mat, I=I, J=J),
                n.chains=2)
mod2_2<-jags.model("M2s2539250.jags", data = list(eff=effect_mat, I=I, J=J), 
                 n.chains=2)

#Use DIC to calculate the deviances of the models.
dic.samples(mod_1,n.iter=30000, n.burnin=7500)
dic.samples(mod2_2,n.iter=50000, n.burnin=20000)

#For the selected run, the penalized deviance of the first model, which has one 
#set of mu values for all 5 of the experiments, is 590.6, while the penalized 
#deviance of the second model, which has a set of mu values for each of the 
#experiments, is 576.4. (The mean deviance of the first is 581.9 with a penalty 
#of 8.71, while the mean deviance of the second is 559.6 with a penalty of 
#16.9). We want a lower penalized deviance for a better model; therefore, the 
#second model should be used. The penalized deviance accounts for the fact that
#we use more parameters in the second model; despite this, the performance is
#better by enough that we should still use this model.

#dev.off()





