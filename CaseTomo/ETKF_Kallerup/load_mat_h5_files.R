#################################################################
### loading results form matlab for Kallerup tomography case ###
# Extracting necessary variables and results from hdf5 files 

# needed library
library(R.matlab)
library(rhdf5)
#options(max.print=1000000)

# File containing both prior, simulated data and rejections sampler posterior realizations
filename_rej <- "rejection/caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G0_rejection_N1000000_di1_out.h5"

# File containing MCMC samples
filename_MH <- "metropolis/caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G0_metropolis_N1000000_di1.mat"


# Data without noise - D
h5_Kallerup_D <- h5read(file = filename_rej, name = "D", index = list(NULL, 1:1000))

# Velocity model - M, prior
h5_Kallerup_M <- h5read(file = filename_rej, name = "M", index = list(NULL, NULL, 1:1000))

# posterior realisations from local and non-local rejection sampler
h5_loc_rej_samples <- h5read(file =filename_rej, name ="post_reals_local")
h5_rej_samples <- h5read(file =filename_rej, name ="post_reals")

# Samples from MH
h5_Kallerup_MCMC <- h5read(file = filename_MH, name = "/reals")

# After all the needed realizations are extracted from the large files
# they are saved into rds objects that can be loaded again into R

saveRDS(h5_Kallerup_D, file= "h5_Kallerup_D_1000.rds")
saveRDS(h5_Kallerup_M, file= "h5_Kallerup_M_1000.rds")


saveRDS(h5_loc_rej_samples, file = "h5_loc_rej_samples.rds") 
saveRDS(h5_rej_samples, file = "h5_rej_samples.rds")
saveRDS(h5_Kallerup_MCMC, file = "h5_Kallerup_MCMC.rds")





