library(R.matlab)
# Installing the rhdf5 package
# install.packages("BiocManager")
# BiocManager::install("rhdf5")
library(rhdf5)
library(fields)
library(viridis)
# Line intersections
library(PlaneGeometry)
library(ggplot2)

source("CaseTomo/functions_support.R")


# Loading the Kallerup data
kallerup <- readMat("CaseTomo/Kallerup/KallerupJensenOutput_small.mat")



kallerup_data_rownames <- rownames(kallerup$data[[1]][[1]])
kallerup_data <- kallerup$data[[1]][[1]]
# different errors
kallerup_obserr_mean <- kallerup$data[[1]][[1]][[7]]
kellerup_obserr_cov <- kallerup$data[[1]][[1]][[8]]
kellerup_CD<- kallerup$data[[1]][[1]][[9]]


# Observed traveltimes
obs_data <- kallerup_data[[4]] 

# Locations of the senders (x,y)
pos_s <- kallerup$ant.pos[,1:2]

# Locations of the receivers (x,y)
pos_r <- kallerup$ant.pos[,3:4]

# getting an overview over content of the hdf5 file
h5ls("CaseTomo/caseTomo_Kallerup_dx25_Feikonal-none_ME0_lu_N200002.h5")

#filename <- "CaseTomo/caseTomo_Kallerup_dx25_Feikonal-none_ME0_lu_N200002.h5"
filename <- "CaseTomo/caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0_lu_noise_N200001.h5"

# Data without noise - D1
h5_Kallerup_D1 <- h5read(file = filename, name = "D1")

# Data with noise - D2
h5_Kallerup_D2 <- h5read(file = filename, name = "D2")

# Velocity model - M1
h5_Kallerup_M1 <- h5read(file = filename, name = "M1")


########### Changable values ############
# Starting and ending values for range of depth and horizontal displacement
# 1:2 horisontal, 3:4 vertical
ax <- c(-0.5,4.0,0,7)
# Spacing between grid nodes
dx <- 0.10
#dx <- 0.25

# Number of ensembles/realisations
n_e <- 100
# Where is the starting point(index) in picking the realisations from the available
r_st <- 160000
# Number of the prior realisation to be plotted
pr_nr <- r_st -1 + 50
# Which posterior realisation should be plotted (should be same as prior)
ps_ne <- 50

# width of the localisation patches
dx_p <- 1.5
dy_p <- 1.75
#dy_p <- 3.5
#########################################



x_coord <- seq(ax[1],ax[2], by = dx)
y_coord <- seq(ax[3],ax[4], by = dx)

nr <- length(x_coord)
nc <- length(y_coord)



# Making of the grid for plotting
grid_list <- as.matrix(expand.grid(x_coord,y_coord))

# Reshaping one of the realisations of the prior model 
reshape_M1_pr <- (matrix(h5_Kallerup_M1[,pr_nr], nrow = nr, ncol = nc,byrow = T))
# reshaping it the second time for plotting 
# mirroring if needed??
#reshape_M1_pr <-  reshape_M1_pr[nrow(reshape_M1_pr):1,]





mean_prior <- matrix((rowMeans(h5_Kallerup_M1)), nrow = nr, ncol = nc,byrow = T)
# Plotting one of the chosen prior realisations
quilt.plot(grid_list, reshape_M1_pr, nx = nr, ny = nc, ylim = c(7,0), ylab = "Y[m]", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(0.07,0.18), asp  = 1)

quilt.plot(grid_list, mean_prior, nx = nr, ny = nc, ylim = c(7,0), ylab = "Y[m]", xlab = "X[m]", legend.lab = "V[m/ns]",  useRaster = T, asp = 1.5)

# Picking out 100 realisations (at random?) to use for the LETKF

x_grid <- matrix(grid_list[,1], nrow = nr, ncol = nc)
y_grid <- matrix(grid_list[,2], nrow = nr, ncol = nc)

image.plot(x_grid, (y_grid), reshape_M1_pr)
quilt.plot(grid_list, reshape_M1_pr, nx = nr, ny = nc,ylab = "Y[m]", xlab = "X[m]", legend.lab = "V[m/ns]")


# "ray tracing": do an apply of lm to the vectors of locations for sources and 
# recievers 

sr_x_coord <- cbind(pos_s[,1], pos_r[,1])
sr_y_coord <- cbind(pos_s[,2], pos_r[,2])

sr_coords <- as.data.frame(cbind(sr_x_coord, sr_y_coord))

colnames(sr_coords) <- c("xS","xR", "yS", "yR")



# intercepts and gradients for each sender-reciever pair
# which also correspond to all of the 412 traces 
#ray_line_eq <- apply(cbind(pos_s,pos_r), 1,find.line)

#quilt.plot(grid_list, reshape_M1_vec_pr, nx = nr, ny = nc, ylim = c(7,0), ylab = "Y[m]", xlab = "X[m]", legend.lab = "V[m/ns]", col = viridis((256)), zlim = c(0.06,0.16))

#segments(pos_r[,1],pos_r[,2],pos_s[,1], pos_s[,2])

##### Identifying which rays cross each of the patches #####
# Logging which trace belongs to which localisation patch #

x_range <- ax[1:2]
y_range <- ax[3:4]


# to make this more robust, do something similar as for alvheim patching 
x_cords <- seq(x_range[1], x_range[2], dx_p)
y_cords <- seq(y_range[1], y_range[2], dy_p)


horis_sides <- as.matrix(expand.grid(x_cords, y_cords))

vert_sides <- as.matrix(expand.grid(y_cords, x_cords))

# for loop for making the matrix of coordinates for horizontal sides
# and matrix of coordinates for vertical sides
# coordinates for start and end of a line giving 4 columns
nLine <- 4
nSides_h <- (length(x_cords)-1)*length(y_cords)

nSides_v <- (length(y_cords)-1)*length(x_cords)

horis_sides_coord <- matrix(NA, nrow = nSides_h, ncol = nLine)
vert_sides_coord <- matrix(NA, nrow = nSides_v, ncol = nLine)

dY <- (length(y_cords)-1)
dX <- (length(x_cords)-1)

for(j in 1:length(y_cords)){
  for(i in 1:(length(x_cords)-1)){
    
    horis_sides_coord[((i)+(j-1)*dX),1:2] <- horis_sides[((i)+(j-1)*(dX+1)),] 
    horis_sides_coord[((i)+(j-1)*dX),3:4] <- horis_sides[((i+1)+(j-1)*(dX+1)),]
    print(((i)+(j-1)*(dX)))
  }
  
}

for(j in 1:length(x_cords)){
  for(i in 1:(length(y_cords)-1)){
    
    vert_sides_coord[((i)+(j-1)*dY),2:1] <- vert_sides[((i)+(j-1)*(dY+1)),] 
    vert_sides_coord[((i)+(j-1)*dY),4:3] <- vert_sides[((i+1)+(j-1)*(dY+1)),]
    print(((i)+(j-1)*dY))
  }
  
}

# List of lists containing indexes of all rays corresponding each of the traveltimes
# crossing different patches
crossing_patch <- list()

patch_nr <- 1

for(j in 1:dY){
  for(i in 1:dX){
    crossing_patch[[patch_nr]] <- list()
    
    # Left side: Line1, Right sides: Line2
    vertLine1 <- Line$new(A = vert_sides_coord[((i-1)*dY+j),1:2], B = vert_sides_coord[((i-1)*dY+j),3:4], FALSE, FALSE)
    vertLine2 <- Line$new(A = vert_sides_coord[((i)*dY+j),1:2],  B = vert_sides_coord[((i)*dY+j),3:4], FALSE, FALSE)
    #print(c(vertLine1,vertLine2))
    # Upper: Line1, Lower: Line2
    horisLine1 <- Line$new(A = horis_sides_coord[((j-1)*dX+i),1:2], B = horis_sides_coord[((j-1)*dX+i),3:4], FALSE, FALSE)
    horisLine2 <- Line$new(A = horis_sides_coord[((j)*dX+i),1:2], B = horis_sides_coord[((j)*dX+i),3:4], FALSE, FALSE)
    #print(c(horisLine1,horisLine2))
    # for loop going checking if any of the traces crosses the line
    # then saves each traces that crosses any of the lines
    # just 3 traces first, change to length(obs_data)
    for(k in 1:length(obs_data)){
      check_cross <- c()
      traceLine <- Line$new(A = pos_s[k,], B = pos_r[k,],FALSE, FALSE)
      
      check_cross <- !is.null(intersectionLineLine(traceLine, vertLine1, strict = TRUE))
      check_cross <- c(check_cross, !is.null(intersectionLineLine(traceLine, vertLine2, strict = TRUE)))
      check_cross <- c(check_cross, !is.null(intersectionLineLine(traceLine, horisLine1, strict = TRUE)))
      check_cross <- c(check_cross, !is.null(intersectionLineLine(traceLine, horisLine2, strict = TRUE)))
      
      if(any(check_cross)==TRUE){
        crossing_patch[[patch_nr]] <- append(crossing_patch[[patch_nr]], k)
      }
      
      
      
    }
    crossing_patch[[patch_nr]]<- unlist(crossing_patch[[patch_nr]])
    patch_nr <- patch_nr + 1
    
  }
}


#### Logging which prior node is belonging to which patch ####
# Extract locations which fall into domain of each of the patches
# For a certain configuration of the patches, this needs to be done only once

# locations of the prior realisations: grid_list

# storage of the indices
prior_patch_indx <- list()

# outer loop, looping over y edge coordinates
patch <- 1
for(i in 1:dY){
  # inner loop, looping over x edge coordinates
  y_edge <- vert_sides_coord[cbind(c(i,i),c(2,4))]
  for(j in 1:dX){
    x_edge <- horis_sides_coord[cbind(c(j,j),c(1,3))]
    prior_patch_indx[[patch]] <- list()
    # inside inner loop have an if-clause checking the coordinates
    
    prior_patch_indx[[patch]] <- which(grid_list[,1]>=x_edge[1] & grid_list[,1]<=x_edge[2] &grid_list[,2]>=y_edge[1] & grid_list[,2]<=y_edge[2])
    
    patch <- patch + 1
    
  }
  
}



######## LETKF step #########
# Prior realisations 

# r_matrix (observation covariance)

# y_matrix - simulated data ( )

Y_avg <- rowMeans(h5_Kallerup_D2[,r_st:(r_st-1+n_e)])

# Anomaly matrix 
Y_matrix <- h5_Kallerup_D2[,r_st:(r_st-1+n_e)] - Y_avg

# Prior: X-matrix

X_avg <- rowMeans(h5_Kallerup_M1[,r_st:(r_st-1+n_e)])
X_matrix <- h5_Kallerup_M1[,r_st:(r_st-1+n_e)] - X_avg

# Temporary R-matrix 
R_temp <- diag(kellerup_CD)
post_velocity <- matrix(NA, nrow = length(X_avg), ncol = n_e)

# Indexes for extraction of variables for each of the patches done in domain_div_loc.R

# Looping through the patches
# need crossing patches list and prior_patch_indx
for(p in 1: length(crossing_patch)){

  obs_data_loc <- obs_data[crossing_patch[[p]],1]
  
  Y_matrix_loc <- Y_matrix[crossing_patch[[p]],]
  
  X_matrix_loc <- X_matrix[prior_patch_indx[[p]],]
  
  y_avg_loc <- Y_avg[crossing_patch[[p]]]
  x_avg_loc <- X_avg[prior_patch_indx[[p]]]
  
  mR_loc <- diag(R_temp[crossing_patch[[p]]],length(crossing_patch[[p]]))
  
  x_final <- loc_entkf(obs_data_loc, mR_loc, Y_matrix_loc, X_matrix_loc, x_avg_loc, y_avg_loc)
  
    
  post_velocity[prior_patch_indx[[p]],] <- x_final
  #browser()
  
}


reshape_M1 <- matrix(post_velocity[,ps_ne], nrow = nr, ncol = nc, byrow=T)
# reshaping it the second time for plotting 
reshape_M1_vec <- matrix(reshape_M1, byrow = T)

# Plotting one of the chosen prior realisations
quilt.plot(grid_list, reshape_M1_vec, nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(0.05,0.17), asp=1)

#image.plot(grid_list, reshape_M1)

reshape_post <- matrix(rowMeans(post_velocity), nrow = nr, ncol = nc, byrow=T)


quilt.plot(grid_list, reshape_post, nx = nr, ny = nc,  ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", legend.lab = "V[m/ns]",zlim = c(0.05, 0.17), asp =1)

#test_data <- expand.grid(y_coord,x_coord)
#test_data$z <- rowMeans(post_velocity)

#ggplot(test_data, aes(Var2, Var1, fill= z)) + 
#  geom_tile()+ylim(7,0)

# test_data_prior <- expand.grid(y_coord,x_coord)
#test_data_prior$z <- h5_Kallerup_M1[,pr_nr]

#ggplot(test_data_prior, aes(Var2, Var1, fill= z)) + 
#  geom_tile()+ylim(7,0)










