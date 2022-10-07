# Distributed ComBat (dComBat)
### Implementation of ComBat for a distributed data setting

Maintained by Andrew Chen, andrewac@pennmedicine.upenn.edu

## 1. Background
ComBat is a widely-used harmonization method that has proven to be effective in both genomic and neuroimaging contexts. For clinical data housed in separate locations however, the ComBat method could not be applied. We adapt ComBat for this distributed data setting as distributed ComBat (dComBat) and provide an implementation based on the current R package for ComBat (https://github.com/Jfortin1/ComBatHarmonization).

If you use this method please cite the following paper:

Chen, A. A., Luo, C., Chen, Y., Shinohara, R. T., & Shou, H. (2022). Privacy-preserving harmonization via distributed ComBat. *NeuroImage*, 248, 118822. https://doi.org/10.1016/j.neuroimage.2021.118822

## 2. Usage
This code is meant to be used without installation to avoid potential complications coordinating across separate data locations. The only current dependency is the package `matrixStats`, but this may be changed in future versions. *neuroCombat_helpers.R* and *neuroCombat.R* are sourced directly from the R implementation of ComBat (https://github.com/Jfortin1/ComBatHarmonization).

Two sample codes are provided *dCombat_central_sample.R* and *dCombat_site_sample.R*. The best way to use our code as follows:

1. Send *dCombat_site_sample.R* to each site and have individual data coordinators adapt that code for their data. The site script will output deidentified summary statistics that can then be sent to a central location.
2. Share the summary statistics with a central location, which modifies *dCombat_central_sample.R* to include these files. This will produce another file, which needs to be sent back to the sites for a second step.
3. After sharing the central location output file, have data coordinators adapt the *dCombat_site_sample.R* code to run a second step, which outputs another set of summary statistics for the last step.
4. Share the second set of summary statistics with the central location, which updates *dCombat_central_sample.R* to output a final set of harmonization parameters.
5. Send these harmonization parameters to each site, which can then perform the final dComBat locally.

## 3. In other programming languages
The Python version of this code is available in the `Python/` directory.
