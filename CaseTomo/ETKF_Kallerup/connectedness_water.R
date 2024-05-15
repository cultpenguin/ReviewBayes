###################################################
########## Integrated values for Kallerup case ####
# Loading the packages
library(fields)
library(paletteer)

# Color palette for plotting
col_p <- paletteer_c("ggthemes::Orange-Blue Diverging", 256, direction = -1)
# Size of the ensemble/number of realizations
n_e <- 500

# Loading the necessary data files
# synthetic data
h5_Kallerup_D <- readRDS(file= "results_all/h5_Kallerup_D_1000.rds")
# prior velocity realizations
h5_Kallerup_M <- readRDS(file= "results_all/h5_Kallerup_M_1000.rds")
# posterior LERS realizations
h5_loc_rej_samples <- readRDS(file = "results_all/h5_loc_rej_samples.rds") 
# posterior rejection realizations
h5_rej_samples <- readRDS(file = "results_all/h5_rej_samples.rds")
# posterior MCMC realizations
h5_Kallerup_MCMC <- readRDS(file = "results_all/h5_Kallerup_MCMC.rds")

# Loading the posterior LETKF realizations
letkf_results <- matrix(unlist(read.table("results_all/final_updated_x_ns_letkf_TEST.txt")),ncol = n_e, nrow = 3266, byrow = F)


##### Setting up the necessary variables for the grid for plotting and reshaping #####

# Extent of the grid
ax <- c(-0.5,4.0,0,7)
# Spacing between grid nodes
dx <- 0.10


x_coord <- seq(ax[1],ax[2], by = dx)
y_coord <- seq(ax[3],ax[4], by = dx)

# dimension of the x-y coordinate system
nr <- length(x_coord)
nc <- length(y_coord)


grid_exp <- expand.grid(y = y_coord,x =x_coord)

xmin <- -2.5
xmax <- 6
ymin <- -2
ymax <- 8

asp_ratio = (xmax-xmin)/(ymax-ymin)


# Function for connectedness
library(raster)
# Function used to identify all the connected areas of the variable of interest which is exceeding a given threshold
connectedness_plot_prep <- function(fluid_transf,inline_small, xline_small, thres, fluid_type, method_type){
  library(RColorBrewer)
  graph_fluid <- matrix(NA, nrow = nrow(fluid_transf), ncol = ncol(fluid_transf))
  for(i in 1:nrow(fluid_transf)){
    for(j in 1:ncol(fluid_transf)){
      if(fluid_transf[i,j]>=thres){
        graph_fluid[i,j] <- 1
      }else{
        graph_fluid[i,j] <- 0
      }
    }
  }
  # Convert the matrix to a raster
  r <- raster(graph_fluid)
  # Identify the clusters
  cl <- clump(r)
  # Get the sizes of the clusters
  clump_sizes <- freq(cl, useNA='no')
  # transposing so that it fits with the original plots
  cl_turned <- (cl)
  
  sorted_clump_sizes <- clump_sizes[order(-clump_sizes[,"count"]), ]
  
  rc_clumped <- cl_turned
  #if(length(unique(clump_sizes[,"count"]))<=11){
  #  col_pal <- brewer.pal(n =11, name = "RdBu")
  #}else{
  #col_pal <- rainbow(length(unique(clump_sizes[,"count"])))
  col_pal <- rainbow(6)
  #}
  
  values(rc_clumped) <- sorted_clump_sizes[match(getValues(cl_turned), sorted_clump_sizes[,"value"]),"count"]
  # Display the labeled raster with the color palette@
  xmin <- -2.5
  xmax <- 6
  ymin <- -2
  ymax <- 8
  
  asp_ratio = (xmax-xmin)/(ymax-ymin)
  # Plotting the connectedness plot
  png(file=paste("Figures/","method", method_type, fluid_type, "_connectedness_plot.png", sep=""),units="px", width = 743*asp_ratio, height  = 743, pointsize = 28)
  plot(rc_clumped, axes = F, col = col_pal, xlab = "X[m]", ylab = "Y[m]", zlim=c(0,750))
  
  axis(1, at=seq(0,1,length.out=5), labels =seq(-0.5,4,length.out=5))
  axis(2,at=seq(0,1,length.out=8),labels = seq(7,0,-1), las = 2)
  
  dev.off()
  # Returning all the sizes of connected areas
  return(clump_sizes)
}




connectedness_area_size <- function(fluid_transf,x_coord, y_coord, thres, indx_small_row=1:which(y_coord==3.0, arr.ind = T)){
  fluid_transf_grid_full <- matrix(fluid_transf, nrow= length(y_coord), ncol = length(x_coord))
  # Specifying the smaller area where we are analyzing the sizes of connectedness
  # For the Kallerup case we divide the area into two
  #ind_small_row <- 1:which(y_coord==3.0, arr.ind = T)
  indx_small_col <- 1:length(x_coord)
  fluid_transf_grid <- fluid_transf_grid_full[indx_small_row, indx_small_col]
  library(RColorBrewer)
  #browser()
  graph_fluid <- matrix(NA, nrow = nrow(fluid_transf_grid), ncol = ncol(fluid_transf_grid))
  for(i in 1:nrow(fluid_transf_grid)){
    for(j in 1:ncol(fluid_transf_grid)){
      if(fluid_transf_grid[i,j]>=thres && !is.nan(fluid_transf_grid[i,j])){
        graph_fluid[i,j] <- 1
      }else{
        graph_fluid[i,j] <- 0
      }
    }
  }
  
  # Convert the matrix to a raster
  r <- raster(graph_fluid)
  # Identify the clusters
  cl <- clump(r)
  # Get the sizes of the clusters
  clump_sizes <- freq(cl, useNA='no')
  # transposing so that it fits with the original plots
  if(nrow(clump_sizes)!=1){
    sorted_clump_sizes <- clump_sizes[order(-clump_sizes[,"count"]), ]
  }else{
    sorted_clump_sizes <- clump_sizes
  }
  
  
  
  return(sorted_clump_sizes[,2])
}


# Threshold of velocity for sand is v>0.1
thres <- 0.1

# Mean of each of the results of the realizations
mean_letkf <- matrix(rowMeans(letkf_results), nrow=71, ncol=46)
std_letkf <- matrix(apply(letkf_results,1,sd), nrow = 71, ncol=46)

mean_mcmc <- rowMeans(h5_Kallerup_MCMC, dims = 2)
std_mcmc <- apply(h5_Kallerup_MCMC, c(1,2), sd)

loc_rej_mean <- rowMeans(h5_loc_rej_samples, dims =2)
loc_rej_std <- apply(h5_loc_rej_samples, c(1,2), sd)

rej_mean <- rowMeans(h5_rej_samples, dims =2)
rej_std <- apply(h5_rej_samples, c(1,2), sd)



##### LERS Accepted realizations #####
# for each posterior realization find connected area
# Above a velocity threshold of 0.1 we know we are dealing with sand
# Calculating and plotting the connectedness plots 
connect_water_LERS <- connectedness_plot_prep(loc_rej_mean, grid_exp$x, grid_exp$y, thres= 0.1, fluid_type="water", method_type="lers")
connect_area_water_LERS_upper <- list()
connect_area_water_LERS_lower <- list()
# for loop for density plot of size
for(i in 1:400){
  # Separate matrices for upper and lower part of the domain
  
  connect_area_water_LERS_upper[[i]] <- connectedness_area_size(h5_loc_rej_samples[,,i],x_coord, y_coord, thres, indx_small_row=1:which(y_coord==3.0, arr.ind = T))
  connect_area_water_LERS_lower[[i]] <- connectedness_area_size(h5_loc_rej_samples[,,i],x_coord, y_coord, thres, indx_small_row=which(y_coord==3.0, arr.ind = T):length(y_coord))
}

# Finding the largest connected area for each section
largest_water_LERS_upper <- unlist(lapply(connect_area_water_LERS_upper, function(x) x[1]), use.names=FALSE)
largest_water_LERS_lower <- unlist(lapply(connect_area_water_LERS_lower, function(x) x[1]), use.names=FALSE)

# Finding the density
den_largest_water_LERS_upper <- density(largest_water_LERS_upper)
den_largest_water_LERS_lower <- density(largest_water_LERS_lower)


#### plain rejection sampler accepted realizations ######
connect_water_rej <- connectedness_plot_prep(rej_mean, grid_exp$x, grid_exp$y, thres= 0.1, fluid_type="water", method_type="rej")
connect_area_water_rej_upper <- list()
connect_area_water_rej_lower <- list()

# For-loop for density plot of size
for(i in 1:15){
  connect_area_water_rej_upper[[i]] <- connectedness_area_size(h5_rej_samples[,,i],x_coord, y_coord, thres, indx_small_row=1:which(y_coord==3.0, arr.ind = T))
  connect_area_water_rej_lower[[i]] <- connectedness_area_size(h5_rej_samples[,,i],x_coord, y_coord, thres, indx_small_row=which(y_coord==3.0, arr.ind = T):length(y_coord))
  
  
}
  

largest_water_rej_upper<- unlist(lapply(connect_area_water_rej_upper, function(x) x[1]), use.names=FALSE)
largest_water_rej_lower <- unlist(lapply(connect_area_water_rej_lower, function(x) x[1]), use.names=FALSE)


den_largest_water_rej_upper <- density(largest_water_rej_upper)
den_largest_water_rej_lower <- density(largest_water_rej_lower)

###### MH Accepted realizaitons ######
connect_water_MH <- connectedness_plot_prep(mean_mcmc, grid_exp$x, grid_exp$y, thres= 0.1, fluid_type="water", method_type="MH")

connect_area_water_MH_upper <- list()
connect_area_water_MH_lower <- list()
# For-loop for density plot of size
for(i in 1:100){
  connect_area_water_MH_upper[[i]] <- connectedness_area_size(h5_Kallerup_MCMC[,,i],x_coord, y_coord, thres, indx_small_row=1:which(y_coord==3.0, arr.ind = T))
  connect_area_water_MH_lower[[i]] <- connectedness_area_size(h5_Kallerup_MCMC[,,i],x_coord, y_coord, thres, indx_small_row=which(y_coord==3.0, arr.ind = T):length(y_coord))

}

largest_water_MH_upper<- unlist(lapply(connect_area_water_MH_upper, function(x) x[1]), use.names=FALSE)
largest_water_MH_lower <- unlist(lapply(connect_area_water_MH_lower, function(x) x[1]), use.names=FALSE)


den_largest_water_MH_upper <- density(largest_water_MH_upper)
den_largest_water_MH_lower <- density(largest_water_MH_lower)


###### LETKF accepted realizations ######
connect_water_LETKF <- connectedness_plot_prep(mean_letkf, grid_exp$x, grid_exp$y, thres= 0.1, fluid_type="water", method_type="letkf")
connect_area_water_LETKF_upper <- list()
connect_area_water_LETKF_lower <- list()
# For-loop for density plot of size 
for(i in 1:500){
  connect_area_water_LETKF_upper[[i]] <- connectedness_area_size(letkf_results[,i],x_coord, y_coord, thres, indx_small_row=1:which(y_coord==3.0, arr.ind = T))
  connect_area_water_LETKF_lower[[i]] <- connectedness_area_size(letkf_results[,i],x_coord, y_coord, thres, indx_small_row=which(y_coord==3.0, arr.ind = T):length(y_coord))
  
  
}

largest_water_LETKF_upper<- unlist(lapply(connect_area_water_LETKF_upper, function(x) x[1]), use.names=FALSE)
largest_water_LETKF_lower <- unlist(lapply(connect_area_water_LETKF_lower, function(x) x[1]), use.names=FALSE)


den_largest_water_LETKF_upper <- density(largest_water_LETKF_upper)
den_largest_water_LETKF_lower <- density(largest_water_LETKF_lower)



#### Upper and lower section plots

png(file=paste("Figures/", "density_connected_sizes_upper_lower_Kallerup.png", sep=""),units="px", width = 2*asp_ratio*2000, height  = 1400, pointsize = 50)

par(mfrow=c(1,2), mar=c(4,4,1,3.5), oma =c(0,0,0,0))

plot(den_largest_water_LERS_upper, main = "", xlab= expression(paste("Number of nodes connected where V>0.1")), lwd =7, xlim = c(0,1100), lty=1, ylim =  c(0,0.01))
lines(den_largest_water_MH_upper, col = "green", lwd = 7, lty=1)
lines(den_largest_water_LETKF_upper, col="purple", lwd= 7, lty=1)
legend(x = "topright",legend = c("LERS", "MH", "LETKF"), col = c("black", "green", "purple"), lty = 1, lwd =7)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -250)


plot(den_largest_water_LERS_lower, main = "", xlab= expression(paste("Number of nodes connected where V>0.1")), lwd =7, xlim = c(0,600), ylim =c(0,0.01))
lines(den_largest_water_MH_lower, col = "green", lwd=7)
lines(den_largest_water_LETKF_lower, col="purple", lwd = 7)
legend(x = "topright",legend = c("LERS", "MH", "LETKF"), col = c("black", "green", "purple"), lty = 1, lwd =7)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -250)


dev.off()






