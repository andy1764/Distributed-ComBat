# Sample code for distributed ComBat, central location

# Make sure to set working directory to source file location!
library(matrixStats)

source("distributedCombat.R")
source("neuroComBat_helpers.R")
source("neuroComBat.R")

# Include outputs from individual sites, can include any number of sites
# Make sure to include site outputs that are on the same step

distributedCombat_central(c("site1_step1.Rdata", "site2_step1.Rdata"), 
                          file = "central_step1.Rdata")

distributedCombat_central(c("site1_step2.Rdata", "site2_step2.Rdata"), 
                          file = "central_step2.Rdata")
