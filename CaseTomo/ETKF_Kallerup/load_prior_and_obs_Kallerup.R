#### Loading:
#### Loading of these files, specificly the prior hdf5 file should be done on a server or from terminal, as
#### RStudio usually doesnt have enough working memory to load such a large file (30gb)
#### - Kallerup prior and corresponding synthetic observations - prior_and_data
#### - Loading the real data (observations) - kallerup
####
library(R.matlab)
library(rhdf5)
# Loading Kallerup data
kallerup <- readMat("Kallerup/caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G0.mat")

# If TRUE - plotting all functions and saving them
#plot_and_save <- FALSE

kallerup_data_rownames <- rownames(kallerup$D)

# Loading kallerup observations
obs_data <- kallerup$D[[3]]
# Loading the observation errors
kallerup_Ct<- kallerup$D[[6]]

# Loading the locations of sources and receivers
kallerup_sources <- kallerup$forward[[5]]
kallerup_receivers <- kallerup$forward[[6]]

# loading the simulated observations
h5_Kallerup_D <- readRDS("results_all/h5_Kallerup_D_1000.rds")
# loading the prior observation/velocity model realisations
h5_Kallerup_M <- readRDS("results_all/h5_Kallerup_M_1000.rds")

