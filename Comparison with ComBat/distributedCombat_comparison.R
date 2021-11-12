library(neuroCombat)
library(matrixStats)
source("../distributedCombat.R")
source("../neuroComBat_helpers.R")
source("../neuroComBat.R")

# Simulate data
set.seed(8888)

p=10000 # Number of features
q=3 # Number of covariates
n=100
batch = rep(1:4, n/4) #Batch variable for the scanner id
batch <- as.factor(batch)

mod = matrix(runif(q*n), n, q) #Random design matrix
dat = matrix(runif(p*n), p, n) #Random data matrix

#### Ground truth ComBat outputs ####
com_out <- neuroCombat(dat, batch, mod)
com_out_ref <- neuroCombat(dat, batch, mod, ref.batch = "1")

#### Distributed ComBat: No reference batch ####
### Step 1
site.outs <- NULL
for (b in unique(batch)) {
  s <- batch == b
  df <- dat[,s]
  bat <- batch[s]
  x <- mod[s,]
  
  site.outs <- c(site.outs, list(distributedCombat_site(df, bat, x)))
}
central <- distributedCombat_central(site.outs)

### Step 2
site.outs <- NULL
for (b in unique(batch)) {
  s <- batch == b
  df <- dat[,s]
  bat <- batch[s]
  x <- mod[s,]
  
  site.outs <- c(site.outs, list(distributedCombat_site(df, bat, x,
                                                        central.out = central)))
}
central <- distributedCombat_central(site.outs)

### Compare distributed vs original
site.outs <- NULL
error <- NULL
perror <- NULL # percent difference
for (b in unique(batch)) {
  s <- batch == b
  df <- dat[,s]
  bat <- batch[s]
  x <- mod[s,]

  site.out <- distributedCombat_site(df, bat, x, central.out = central)
  site.outs <- c(site.outs, site.out)
  
  error <- c(error, max(c(com_out$dat.combat[,s] - site.out$dat.combat)))
  perror <- c(perror,
              max(c(abs(com_out$dat.combat[,s] - site.out$dat.combat)/
                      site.out$dat.combat)))
}

error # Maximum difference by batch
perror # Maximum difference as percent by batch

#### Distributed ComBat: With reference batch ####
### Step 1
site.outs <- NULL
for (b in unique(batch)) {
  s <- batch == b
  df <- dat[,s]
  bat <- batch[s]
  x <- mod[s,]
  
  site.outs <- c(site.outs, list(distributedCombat_site(df, bat, x, ref.batch = "1")))
}
central <- distributedCombat_central(site.outs, ref.batch = "1")

### Step 2
site.outs <- NULL
for (b in unique(batch)) {
  s <- batch == b
  df <- dat[,s]
  bat <- batch[s]
  x <- mod[s,]
  
  site.outs <- c(site.outs, list(distributedCombat_site(df, bat, x, ref.batch = "1",
                                                        central.out = central)))
}
central <- distributedCombat_central(site.outs, ref.batch = "1")

### Compare distributed vs original
site.outs <- NULL
error_ref <- NULL
perror_ref <- NULL # percent difference
for (b in unique(batch)) {
  s <- batch == b
  df <- dat[,s]
  bat <- batch[s]
  x <- mod[s,]
  
  site.out <- distributedCombat_site(df, bat, x, ref.batch = "1", 
                                     central.out = central)
  site.outs <- c(site.outs, site.out)
  
  error_ref <- c(error_ref, max(c(com_out_ref$dat.combat[,s] - site.out$dat.combat)))
  perror_ref <- c(perror_ref,
              max(c(abs(com_out_ref$dat.combat[,s] - site.out$dat.combat)/
                      site.out$dat.combat)))
}

#### Display errors ####
error # Maximum difference by batch
perror # Maximum difference as percent by batch

error_ref # Maximum difference by batch
perror_ref # Maximum difference as percent by batch
