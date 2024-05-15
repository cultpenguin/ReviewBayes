####################################################
##### Independent script for running NS-LETKF ######
####################################################
# Only for the local run
#if(getwd() != "/Users/minasp/Library/CloudStorage/OneDrive-NTNU/Git_projects/ReviewBayes/CaseTomo"){
#  setwd("/Users/minasp/Library/CloudStorage/OneDrive-NTNU/Git_projects/ReviewBayes/CaseTomo")
#}
source("load_prior_and_obs_Kallerup.R")
source("functions_support.R")


# The code is for transformation of the original bimodal variables to
# a unimodal standard normal variable 
set.seed(357)
# Specify the ensemble size
n_e <- 500

# Localization variables 
# width of the localisation patches
dx_p <- 1.5

# depth/length of the localization patches
dy_p <- 1.75

# Where is the starting point(index) in picking the realisations from the available
# prior realizations where #prior realization >= n_e
r_st <- 1

# Locations of the senders (x,y)
pos_s <- kallerup_sources

# Locations of the receivers (x,y)
pos_r <- kallerup_receivers



# Domain and x- and y-coordinates
ax <- c(-0.5,4.0,0,7)
# Spacing between grid nodes
dx <- 0.10


# Making the x- and y- coordinates before the grid is made
x_coord <- seq(ax[1],ax[2], by = dx)
y_coord <- seq(ax[3],ax[4], by = dx)


# grid variable, also used later for plotting
grid_exp <- expand.grid(y = y_coord,x =x_coord)

# alternative way for constructing the grid
grid_alt <- cbind((grid_exp$x),(grid_exp$y))
grid_alt <- as.matrix(expand.grid(y_coord, x_coord))
grid_alt[,1:2] <- grid_alt[,2:1]


colnames(grid_alt) <- c("x","y")

# size of the x and y dimension
nr <- length(x_coord)
nc <- length(y_coord)


# loading the prior realizations from the file
sample_x <- h5_Kallerup_M[,,1:n_e]

# converting the realizations to the dimension that fits with the rest of the code
# from 3D array to 2D array
mid_sample_x <- matrix(sample_x, nrow = prod(dim(sample_x)[1:2]), ncol = dim(sample_x)[3])
# turn into a vector, columns stacked
long_sample_x <- c(mid_sample_x)
# compute empirical cdf
quantiles_sample_x_func <- ecdf(sample_x)
# returns percentiles for the sample_x
quantiles_sample_x <- quantiles_sample_x_func(long_sample_x)
look_up <- cbind(long_sample_x, quantiles_sample_x)

# Returns the z-values from N(0,1)
quantiles_sample_x[which(quantiles_sample_x == 1)] <- 0.9999999999999999
# Converting percentiles to z-scores
z_from_quantiles <- qnorm(quantiles_sample_x)
z_matrix <- matrix(z_from_quantiles, nrow = (nr*nc), ncol = n_e, byrow = F)

# Preparing the transformed variables for the LETKF
Z_avg <- rowMeans(z_matrix)
Z_matrix <- z_matrix - Z_avg

post_velocity_z <- matrix(NA, nrow = length(Z_avg), ncol = n_e)
# Finding the lines that are crossing the localization patches defined for this case
source("line_crossing.R")

# Need the Y-matrix for LETKF 
Y_avg <- rowMeans(h5_Kallerup_D[,r_st:(r_st-1+n_e)])

# Anomaly matrix 
Y_matrix <- h5_Kallerup_D[,r_st:(r_st-1+n_e)] - Y_avg

# Temporary R-matrix 
R_temp <- kallerup_Ct

for(p in 1: length(crossing_patch)){
  
  # extracting the data for the given patch
  obs_data_loc <- obs_data[crossing_patch[[p]],1]
  
  # extracting simulated data for the given patch
  Y_matrix_loc <- Y_matrix[crossing_patch[[p]],]
  
  # extracting prior velocity realizations for the given patch
  Z_matrix_loc <- Z_matrix[prior_patch_indx[[p]],]
  
  y_avg_loc <- Y_avg[crossing_patch[[p]]]
  z_avg_loc <- Z_avg[prior_patch_indx[[p]]]
  
  # Observation error correpsonding to the given patch
  mR_loc <- R_temp[crossing_patch[[p]],crossing_patch[[p]]]
  
  # Final update of the prior realizations in the patch p using LETKF
  z_final <- loc_entkf(obs_data_loc, mR_loc, Y_matrix_loc, Z_matrix_loc, z_avg_loc, y_avg_loc, n_e)
  
  # inputting back the updated patch
  post_velocity_z[prior_patch_indx[[p]],] <- z_final
  #browser()
  
}
# The transformation back from the standard normal unimodal variabe to the bimodal variable

# the pnorm of the variable - we get quantiles there
inv_quantiles <- pnorm(post_velocity_z)

# reshaping
test_quant_inv <- matrix(unname(quantile(sample_x, probs = inv_quantiles, type = 3)),nrow = 3266, ncol = n_e, byrow = T)

# doing some checks and pitfalls
if(any(is.na(inv_quantiles))){warning("There is possibly of a quantile that is 1.")}


tol <- 1e-04
x_updated <- matrix(NA, nrow = length(Z_avg), ncol = n_e)

which_comparison <- function(x,y){}

indices_mx <-  which(round(quantiles_sample_x, digits = 5) == round(inv_quantiles[1,1],digits = 5))
I <- nrow(x_updated)
# matching the quantiles back to the bimodal distribution through ivnerse normal score
start_time <- Sys.time()
for(j in 1:ncol(x_updated)){
  for(i in 1:I){
    
    x_updated[i,j] <-  mean(long_sample_x[which(abs(quantiles_sample_x - inv_quantiles[i,j])<tol)])
    if(is.nan(x_updated[i,j]))
    {print(min(abs(quantiles_sample_x - inv_quantiles[i,j])))}
    
  }
  
}
end_time <- Sys.time()
tot_time <- end_time - start_time
tot_time
# Saving the total time it takes
saveRDS(tot_time, file="tot_time_ns_letkf.rds")

# writing to file the final results
write.table(x_updated, file = "final_updated_x_ns_letkf_TEST.txt", sep = " ", col.names = FALSE, row.names = FALSE)




