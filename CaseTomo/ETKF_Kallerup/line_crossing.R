########################################################
#### Code calculating which lines cross the patches ####
#### in the localization procedure for the Kallerup ####
library(PlaneGeometry)
# Necessary variables for this code: 

# coordinates for start and end of a line giving 4 columns
x_range <- ax[1:2]
y_range <- ax[3:4]



x_cords <- seq(x_range[1], x_range[2], dx_p)
y_cords <- seq(y_range[1], y_range[2], dy_p)

# Defining vertical and horizontal sides of the localization patch 
nLine <- 4
nSides_h <- (length(x_cords)-1)*length(y_cords)

nSides_v <- (length(y_cords)-1)*length(x_cords)

horis_sides <- as.matrix(expand.grid(x_cords, y_cords))

vert_sides <- as.matrix(expand.grid(y_cords, x_cords))


# end point coordinates of horisontal and vertical edges of the patches
horis_sides_coord <- matrix(NA, nrow = nSides_h, ncol = nLine)
vert_sides_coord <- matrix(NA, nrow = nSides_v, ncol = nLine)

dY <- (length(y_cords)-1)
dX <- (length(x_cords)-1)
# Identifying corners in coordinates of each of the localization patches 
for(j in 1:length(y_cords)){
  for(i in 1:(length(x_cords)-1)){
    
    horis_sides_coord[((i)+(j-1)*dX),1:2] <- horis_sides[((i)+(j-1)*(dX+1)),] 
    horis_sides_coord[((i)+(j-1)*dX),3:4] <- horis_sides[((i+1)+(j-1)*(dX+1)),]
    #print(((i)+(j-1)*(dX)))
  }
  
}

for(j in 1:length(x_cords)){
  for(i in 1:(length(y_cords)-1)){
    
    vert_sides_coord[((i)+(j-1)*dY),2:1] <- vert_sides[((i)+(j-1)*(dY+1)),] 
    vert_sides_coord[((i)+(j-1)*dY),4:3] <- vert_sides[((i+1)+(j-1)*(dY+1)),]
    #print(((i)+(j-1)*dY))
  }
  
}

# List of lists containing indexes of all rays corresponding each of the traveltimes
# crossing different patches
crossing_patch <- list()

patch_nr <- 1

for(j in 1:dY){
  for(i in 1:dX){
    crossing_patch[[patch_nr]] <- list()
    #browser()
    
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
    #browser()
    x_edge <- horis_sides_coord[cbind(c(j,j),c(1,3))]
    prior_patch_indx[[patch]] <- list()
    
    # rounding function is added due to numerical comparison of floating point values
    prior_patch_indx[[patch]] <- which(round(grid_alt[,1], digits = 2)>=round(x_edge[1], digits =2) & round(grid_alt[,1], digits = 2)<=round(x_edge[2], digits =2) & round(grid_alt[,2], digits = 2)>=round(y_edge[1], digits = 2) & round(grid_alt[,2], digits = 2)<=round(y_edge[2], digits=2))
    
    
    patch <- patch + 1
    
  }
  
}
