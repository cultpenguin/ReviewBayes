#################################################
#### Kallerup loading and plotting the data #####
#################################################
library(rhdf5)
library(fields)
library(viridis)
library(paletteer)

# Specifying the color palette
col_p <- paletteer_c("ggthemes::Orange-Blue Diverging", 256, direction = -1)
# ensemble size/number of realizations
n_e <- 500

# Loading the synthetic data
h5_Kallerup_D <- readRDS(file= "results_all/h5_Kallerup_D_1000.rds")
# Loading the prior realizations of the velocity model
h5_Kallerup_M <- readRDS(file= "results_all/h5_Kallerup_M_1000.rds")
# Loading the posterior realizations of the LERS
h5_loc_rej_samples <- readRDS(file = "results_all/h5_loc_rej_samples.rds") 
# Loading the posterior realizations of the rejection sampler
h5_rej_samples <- readRDS(file = "results_all/h5_rej_samples.rds")

h5_Kallerup_MCMC <- readRDS(file = "results_all/h5_Kallerup_MCMC.rds")

# NS_LETKF results 
letkf_results <- matrix(unlist(read.table("results_all/final_updated_x_ns_letkf_TEST.txt")),ncol = n_e, nrow = 3266, byrow = F)


### Setting up necessary variables for defining the grid which is needed for plotting and reshaping ###
# Limits of the x-y grid
ax <- c(-0.5,4.0,0,7)
# Spacing between grid nodes
dx <- 0.10

# Specifying the x and y coordinates
x_coord <- seq(ax[1],ax[2], by = dx)
y_coord <- seq(ax[3],ax[4], by = dx)

# dimension of y and x axis
nr <- length(x_coord)
nc <- length(y_coord)

# setting up the grid
grid_exp <- expand.grid(y = y_coord,x =x_coord)
# alternative setup of grid
grid_alt <- cbind((grid_exp$x),(grid_exp$y))
grid_alt <- as.matrix(expand.grid(y_coord, x_coord))
grid_alt[,1:2] <- grid_alt[,2:1]
#########

# Computing the mean and standard deviation of LETKF ensemble
mean_letkf <- rowMeans(letkf_results)
std_letkf <- apply(letkf_results,1,sd)
# Computing the mean and standard deviation of posterior MCMC realizations
mean_mcmc <- rowMeans(h5_Kallerup_MCMC, dims = 2)
std_mcmc <- apply(h5_Kallerup_MCMC, c(1,2), sd)
# Computing the mean and standard deviation of posterior realizaitons LERS
loc_rej_mean <- rowMeans(h5_loc_rej_samples, dims =2)
loc_rej_std <- apply(h5_loc_rej_samples, c(1,2), sd)
# Computing the mean and standard deviation of posterior realizaitons rejection sampler
rej_mean <- rowMeans(h5_rej_samples, dims =2)
rej_std <- apply(h5_rej_samples, c(1,2), sd)

# Limits of z-axis for plotting
min_z_samp <- 0.07
max_z_samp <- 0.25

min_z_mean <- 0.07
max_z_mean <- 0.20

# Calculations of asp. ratio factor for figure size
xmin <- -2.5
xmax <- 6
ymin <- -2
ymax <- 8

asp_ratio = (xmax-xmin)/(ymax-ymin)

#### Posterior plot ####
png(file=paste("Figures/", "posterior_LETKF.png", sep=""),units="px", width = 1400*2*asp_ratio, height  = 1400, pointsize = 50)
par(mfrow=c(1,2), mar=c(4,4,1,4.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, mean_letkf, nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =6.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, std_letkf, nx = nr, ny = nc, xlim =c(-0.5,4), ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", axes = FALSE, col = col_p, legend.shrink = 0.9, legend.line=3, legend.mar = 6.1)
axis(1)
box()
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
dev.off()

#### Plotting the density change ####
indx <- which(grid_exp$x == "2" & grid_exp$y == "2")
indx <- which(grid_exp$x == "0.5" & grid_exp$y == "3")
indx <- which(grid_exp$x == "2" & grid_exp$y == "1.3")

# Coordinates for accessing values from the prior model, due to the way array is constructed
indx_x <- which(x_coord == 2.0)
indx_y <- which(y_coord == 1.3)

den_prio <- density(h5_Kallerup_M[indx_y,indx_x,1:n_e])
den_post <- density(letkf_results[indx,])

png(file=paste("Figures/", "density_change_LETKF.png", sep=""),units="px", width = 1600, height  = 1600*asp_ratio, pointsize = 48)
par(mfrow=c(1,1), mar=c(4,4,3,2))
plot(den_prio, xlim = c(min(den_prio$x,den_post$x), max(den_prio$x,den_post$x)), main = "", xlab ="V[m/ns]", ylab = "Density", col = "black", lwd = 4)
points(x = h5_Kallerup_M[indx_y, indx_x,1:10], y = rep(0.1,10), lwd =4)
lines(den_post, col = "red", lwd =4)
points(x = letkf_results[indx,1:10], y = rep(5,10), col = "red", lwd =4)
arrows(x0 = h5_Kallerup_M[indx_y,indx_x,1:10], y0 = rep(0.1,10), x1 = letkf_results[indx,1:10], y1 = rep(5,10), length=0.06, lwd =2.5)
legend(x = "topright",legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, lwd =3)
dev.off()


## Prior realizations ##
png(file=paste("Figures/", "prior_realizations.png", sep=""),units="px", width = 4*1400*asp_ratio, height  = 1400, pointsize = 100)

par(mfrow=c(1,4), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_M[,,400], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =7.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_M[,,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.9, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_M[,,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar = 7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_M[,,310], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar = 7.1, legend.cex = 0.7)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
dev.off()


#### EXTENDED MH RESULTS ####
png(file=paste("Figures/", "posterior_MCMC.png", sep=""),units="px", width = 2*asp_ratio*1400, height  = 1400, pointsize = 50)
par(mfrow=c(1,2), mar=c(4,4,1,4.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, mean_mcmc, nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", legend.lab = "",  zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =6.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, std_mcmc, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", main = "", col = col_p, legend.shrink = 0.9, legend.line=3, axes =F)
axis(1)
box()
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
dev.off()

#### EXTENDED REJECTION RESULTS ####


png(file=paste("Figures/", "posterior_rej.png", sep=""),units="px", width = 2*asp_ratio*1400, height  = 1400, pointsize = 50)
par(mfrow=c(1,2), mar=c(4,4,1,4.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, rej_mean, nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", legend.lab = "",  zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =6.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, rej_std, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", main = "", col = col_p, legend.shrink = 0.9, legend.line=3, axes =F)
axis(1)
box()
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
dev.off()


png(file=paste("Figures/", "posterior_rej_local.png", sep=""),units="px", width = 2*asp_ratio*1400, height  = 1400, pointsize = 50)
par(mfrow=c(1,2), mar=c(4,4,1,4.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, loc_rej_mean, nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", legend.lab = "",  zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =6.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, loc_rej_std, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", main = "", col = col_p, legend.shrink = 0.9, legend.line=3, axes =F)
axis(1)
box()
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 1, at = -2)
dev.off()

#### Plotting 4 posterior samples of each method ####

# LETKF
png(file=paste("Figures/", "post_realizations_letkf.png", sep=""),units="px", width = 4*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,4), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,400], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =7.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.9, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,310], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3,legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)

dev.off()

# Local rejection sampler

png(file=paste("Figures/", "post_realizations_loc_rej.png", sep=""),units="px", width = 4*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,4), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,400], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =7.1, legend.cex =0.7)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.9, legend.line = 3, legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3,legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,310], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3,legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
dev.off()

# Rejection sampler

png(file=paste("Figures/", "post_realizations_rej.png", sep=""),units="px", width = 4*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,4), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,2], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =7.1, legend.cex =0.7)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,5], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.9, legend.line = 3, legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,11], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,14], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
dev.off()

# MH samples

png(file=paste("Figures/", "post_realizations_mcmc.png", sep=""),units="px", width = 4*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,4), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_MCMC[,,41], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.9, legend.line = 3, legend.mar =7.1, legend.cex =0.7)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_MCMC[,,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.9, legend.line = 3, legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_MCMC[,,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_MCMC[,,77], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.9, legend.line=3, legend.mar =7.1, legend.cex =0.7)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -3)
dev.off()

##############################################################
###### Figures with realizations and mean and std in one ##### 
##############################################################


png(file=paste("Figures/", "post_realizations_mean_std_letkf.png", sep=""),units="px", width = 6*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,6), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,400], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2.5)
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, letkf_results[,310], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3,legend.mar =7.1)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, mean_letkf, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), axes = F, col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("e)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, std_letkf, nx = nr, ny = nc, xlim =c(-0.5,4), ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(0,0.044),axes = FALSE, col = col_p, legend.shrink = 0.8, legend.line=3, legend.mar = 10.1, legend.cex =0.7)
axis(1)
box()
mtext(~ bold("f)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)

dev.off()

# Local rejection sampler
png(file=paste("Figures/", "post_realizations_mean_std_loc_rej.png", sep=""),units="px", width = 6*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,6), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,400], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2.5)
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y,  h5_loc_rej_samples[,,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, h5_loc_rej_samples[,,310], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3,legend.mar =7.1)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, loc_rej_mean, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), axes = F, col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("e)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, loc_rej_std, nx = nr, ny = nc, xlim =c(-0.5,4), ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(0,0.044),axes = FALSE, col = col_p, legend.shrink = 0.8, legend.line=3, legend.mar = 10.1, legend.cex =0.7)
axis(1)
box()
mtext(~ bold("f)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)

dev.off()

# Rejection sampler
png(file=paste("Figures/", "post_realizations_mean_std_rej.png", sep=""),units="px", width = 6*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,6), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,2], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2.5)
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,5], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y,  h5_rej_samples[,,11], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, h5_rej_samples[,,14], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3,legend.mar =7.1)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, rej_mean, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), axes = F, col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("e)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, rej_std, nx = nr, ny = nc, xlim =c(-0.5,4), ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(0,0.044),axes = FALSE, col = col_p, legend.shrink = 0.8, legend.line=3, legend.mar = 10.1, legend.cex =0.7)
axis(1)
box()
mtext(~ bold("f)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)

dev.off()



# MH samples
png(file=paste("Figures/", "post_realizations_mean_std_mcmc.png", sep=""),units="px", width = 6*asp_ratio*1400, height  = 1400, pointsize = 100)

par(mfrow=c(1,6), mar=c(4,4,1,1.5), oma =c(0,0,0,2))
quilt.plot(grid_exp$x, grid_exp$y,h5_Kallerup_MCMC[,,41], nx = nr, ny = nc, ylim = c(7,0),ylab = "Y[m]", xlab = "X[m]", zlim = c(min_z,max_z), main = "", col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
mtext(~ bold("a)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2.5)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_MCMC[,,25], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]",zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("b)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y,   h5_Kallerup_MCMC[,,30], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3, legend.mar =7.1)
axis(1)
mtext(~ bold("c)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, h5_Kallerup_MCMC[,,77], nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), main = "", col = col_p, axes =F, legend.shrink = 0.8, legend.line=3,legend.mar =7.1)
axis(1)
mtext(~ bold("d)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, mean_mcmc, nx = nr, ny = nc, ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "", zlim = c(min_z,max_z), axes = F, col = col_p, legend.shrink=0.8, legend.line = 3, legend.mar =7.1)
axis(1)
mtext(~ bold("e)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)
quilt.plot(grid_exp$x, grid_exp$y, std_mcmc, nx = nr, ny = nc, xlim =c(-0.5,4), ylim = c(7,0),ylab = "", xlab = "X[m]", legend.lab = "V[m/ns]", zlim = c(0,0.044),axes = FALSE, col = col_p, legend.shrink = 0.8, legend.line=3, legend.mar = 10.1, legend.cex =0.7)
axis(1)
box()
mtext(~ bold("f)"), side = 3, line = 0.5, padj = 1, adj = 0, cex = 0.9, at = -2)

dev.off()






