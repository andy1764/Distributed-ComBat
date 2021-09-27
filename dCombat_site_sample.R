# Sample code for distributed ComBat

# Make sure to set working directory to source file location!
# install.packages("matrixStats")
library(matrixStats)

source("distributedCombat.R")
source("neuroComBat_helpers.R")
source("neuroComBat.R")

# You will need the following variables:
#  - dat: features x subject data matrix for this site
#  - bat: batch identifiers, needs to have same factor levels across sites
#  - mod: covariates to protect in the data, usually output of stats:model.matrix

# first, get summary statistics needed for LS estimation
distributedCombat_site(dat, bat, mod, file = "site1_step1.Rdata")

# after step 1 at central site, get summary statistics for sigma estimation
distributedCombat_site(dat, bat, mod, file = "site1_step2.Rdata",
                       central.out = "central_step1.Rdata")

# after step 2 at central site, get harmonized data
distributedCombat_site(dat, bat, mod, file = "site1_harmonized_data.Rdata",
                       central.out = "central_step2.Rdata")
