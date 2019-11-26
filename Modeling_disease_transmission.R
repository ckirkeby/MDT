

# Code for showing simple examples of the different types of models created in R.


#### Difference equation model example ####

beta <- 0.5 # transmission rate
timesteps <- 1:10 # Time steps
S <- 9
I <- 1
N <- S + I
n.inf <- numeric()
n.sus <- numeric()
n.pop <- numeric()

for (i in timesteps) 
{
  dI <- beta * S * I / N
  I <- I + dI
  S <- S - dI
  N <- S + I
  n.inf[i] <- I
  n.sus[i] <- S
  n.pop[i] <- I + S
}

plot(timesteps, n.pop, type="b", ylab="Proportion of population", xlab="Time steps", main="Difference equation model example", ylim=c(0,10))
points(n.inf, col="red", type="b")
points(n.sus, col="blue", type="b")
legend("right", bty="n", lty=b, legend=c("Total population", "No. infected", "No. susceptible"), 
       text.col = c("black", "red", "blue"), col = c("black", "red", "blue"))

#### Differential equation model example ####


## Load deSolve package
library(deSolve)

## Create an SIR function
SI <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I 
    dI <-  beta * S * I
    
    return(list(c(dS, dI)))
  })
}

### Set parameters
# Proportion infected at start:
inf.start <- 0.1
## Proportion in each compartment: 10% infected and 90% susceptible
init       <- c(S = 1-inf.start, I = inf.start)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.5)
## Time frame
timesteps <- c(1:10)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = timesteps, func = SI, parms = parameters)
## change to data frame
out <- as.data.frame(out)
out$N <- out$S + out$I

plot(timesteps, out$N, type="b", xlab="Time steps", 
     main="Differential equation model example", ylim=c(0,1), 
     ylab="Proportion of population")
points(out$I, col="red", type="b")
points(out$S, col="blue", type="b")
legend("right", bty="n", lty=b, legend=c("Total population", "No. infected", "No. susceptible"), 
       text.col = c("black", "red", "blue"), col = c("black", "red", "blue"))


#### Stochastic population-based mechanistic model example ####

beta <- 0.5 # transmission rate
timesteps <- 1:10 # Time steps
S <- 9
I <- 1
N <- S + I
n.inf <- numeric()
n.sus <- numeric()
n.pop <- numeric()

for(i in timesteps)
{
  PI <- 1 - exp( - beta * I / N)
  dI <- sum( rbinom(S , 1 , prob=PI) )
  I <- I + dI
  S <- S - dI
  N <- S + I
  n.inf[i] <- I
  n.sus[i] <- S
  n.pop[i] <- I + S
}

plot(timesteps, n.pop, type="b", ylab="Proportion of population", xlab="Time steps", 
     main="Stochastic population-based mechanistic model example", ylim=c(0,10))
points(n.inf, col="red", type="b")
points(n.sus, col="blue", type="b")
legend("right", bty="n", lty=b, legend=c("Total population", "No. infected", "No. susceptible"), 
       text.col = c("black", "red", "blue"), col = c("black", "red", "blue"))


#### Stochastic individual-based mechanistic model example ####

beta <- 0.5 # transmission rate
timesteps <- 1:10 # Time steps
pop <- data.frame(ID = c(1:10), inf.status = c(1 , rep(0, 9)) )
n.inf <- numeric()
n.sus <- numeric()
n.pop <- numeric()

for(i in timesteps)
{
  I <- length(pop$inf.status[pop$inf.status == 1])
  S <- length(pop$inf.status[pop$inf.status == 0])
  N <- length(pop$inf.status)
  PI <- 1 - exp( - beta * I / N)
  new.inf <- rbinom(S , 1 , prob=PI) 
  pop$inf.status[pop$inf.status == 0] <- new.inf
  n.inf[i] <- I
  n.sus[i] <- S
  n.pop[i] <- I + S
}

plot(timesteps, n.pop, type="b", ylab="Proportion of population", xlab="Time steps", 
     main="Stochastic individual-based mechanistic model example", ylim=c(0,10))
points(n.inf, col="red", type="b")
points(n.sus, col="blue", type="b")
legend("left", bty="n", lty=b, legend=c("Total population", "No. infected", "No. susceptible"), 
       text.col = c("black", "red", "blue"), col = c("black", "red", "blue"))


#### Stochastic individual-based mechanistic simulation model with multiple iterations example ####

MaxIterations <- 10
Output <- matrix(numeric(0),ncol=3)
timesteps <- 1:10 # Time steps

for(j in 1:MaxIterations)
{
  beta <- runif(1,0.4,0.6) # transmission rate
  pop <- data.frame(ID = c(1:10), inf.status = c(1 , rep(0, 9)) )
  n.inf <- numeric()
  n.sus <- numeric()
  n.pop <- numeric()
  
  for(i in timesteps)
  {
    I <- length(pop$inf.status[pop$inf.status == 1])
    S <- length(pop$inf.status[pop$inf.status == 0])
    N <- length(pop$inf.status)
    PI <- 1 - exp( - beta * I / N)
    new.inf <- rbinom(S , 1 , prob=PI) 
    pop$inf.status[pop$inf.status == 0] <- new.inf
    n.inf[i] <- I
    n.sus[i] <- S
    n.pop[i] <- I + S
    
    Output <- rbind(Output,c(i,sum(pop$inf.status==0),sum(pop$inf.status==1)))
  }
}

## If you run the model for a single iteration, you can observe the progress of the infection for that iteration using this code
plot(timesteps, n.pop, type="b", ylab="Y", xlab="Time steps", 
     main="Stochastic individual-based mechanistic model example", ylim=c(0,10))
points(n.inf, col="red", type="b")
points(n.sus, col="blue", type="b")
legend("left", bty="n", lty=b, legend=c("Total population", "No. infected", "No. susceptible"), 
       text.col = c("black", "red", "blue"), col = c("black", "red", "blue"))

## If you run the model for > 1 iteration, then you can observe the effect of randomness on infection using the following code
plot(Output[,1],Output[,2],xlab="Time steps", ylab="Number of animals", ylim=c(0,10), typ="l", col="blue")
lines(Output[,1],Output[,3],col="red")
legend("left", bty="n", lty=b, legend=c("No. infected", "No. susceptible"), 
       text.col = c("red", "blue"), col = c("red", "blue"))

## To observe the progress of infection based on all iterations, the median number can be ploted over time as follows
Sus    <- sapply(unique(Output[,1]),function(x) median(Output[Output[,1]==x,2]))
Infect <- sapply(unique(Output[,1]),function(x) median(Output[Output[,1]==x,3]))
plot(unique(Output[,1]),Sus,xlab="Time steps", ylab="Number of animals",ylim=c(0,10),typ="l", col="blue")
lines(unique(Output[,1]),Infect,col="red")
legend("left", bty="n", lty=b, legend=c("No. infected", "No. susceptible"), 
       text.col = c("red", "blue"), col = c("red", "blue"))


#### Sensitivity analysis example ####

# Make the model (except the beta definition) as a function of beta:
# Here we use the individual-based mechanistic model, but it could be any of those described above.
model <- function(beta)
{
  pop <- data.frame(ID = c(1:10), inf.status = c(1 , rep(0, 9)) )
  n.inf <- numeric()
  n.sus <- numeric()
  n.pop <- numeric()
  
  for(i in timesteps)
  {
    I <- length(pop$inf.status[pop$inf.status == 1])
    S <- length(pop$inf.status[pop$inf.status == 0])
    N <- length(pop$inf.status)
    PI <- 1 - exp( - beta * I / N)
    new.inf <- rbinom(S , 1 , prob=PI) 
    pop$inf.status[pop$inf.status == 0] <- new.inf
    n.inf[i] <- I
    n.sus[i] <- S
    n.pop[i] <- I + S
  }
  
  return(data.frame(n.pop=n.pop, n.inf=n.inf, n.sus=n.sus))  
}

beta <- 0.5 # transmission rate
timesteps <- 1:10 # Time steps

tmp <- model(beta)

plot(timesteps, tmp$n.pop, type="b", ylab="Number of animals", xlab="Time steps", 
     main="Stochastic individual-based mechanistic model example", ylim=c(0,10))
points(tmp$n.inf, col="red", type="b")
points(tmp$n.sus, col="blue", type="b")
legend("left", bty="n", lty=b, legend=c("Total population", "No. infected", "No. susceptible"), 
       text.col = c("black", "red", "blue"), col = c("black", "red", "blue"))

# Now the model can be run several times with a new result because it is stochastic.
# We can then use the model to find the distribution of the number of infected individuals:

# We loop over a number of iterations:
iterations <- 1:1000
infected.end <- numeric()
for(j in iterations)
{
  tmp <- model(beta)
  infected.end[j] <- tmp$n.inf[10]
}
hist(infected.end, xlab="Number of infected individuals at the end of the simulation", main=" ")

# In this case, most of the population is most frequently infected.
table(infected.end)
# We can then decrease the beta and look again:

beta <- 0.35 # transmission rate
# We loop over a number of iterations:
iterations <- 1:1000
infected.end <- numeric()
for(j in iterations)
{
  tmp <- model(beta)
  infected.end[j] <- tmp$n.inf[10]
}
hist(infected.end, main=" ", xlab="Number of infected individuals at the end of the simulation")

# Now the most frequent number of infected animals are around 6. Notice that there
# are a large proportion of iterations where the only infected individual is the 
# one that was infected from the beginning: the epidemic never took off.
# We can see this clearly in the table:
table(infected.end)

#### Convergence ####


beta <- 0.5 # transmission rate
timesteps <- 1:10 # Time steps

# We use the model defined above, in the loop.
iterations <- 1:1000
infected.end <- numeric()
for(j in iterations)
{
  tmp <- model(beta)
  infected.end[j] <- tmp$n.inf[10]
}

# Now find the variance for 1 to 1000 iterations:
conv <- numeric()
for(u in 2:length(iterations))
{
  conv[u] <- var(infected.end[1:u])
}

plot(conv, xlab="Number of iterations", ylab="Variance", type="l")

# Generally, it looks like 400 iterations are sufficient for convergence with the default values used here




# Copyright (c) Carsten Kirkeby 2019
