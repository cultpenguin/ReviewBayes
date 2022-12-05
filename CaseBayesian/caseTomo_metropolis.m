clear all;close all
load caseBayesian_dx15_Fray_2d-none_ME0.mat
load caseBayesian_dx15_Feikonal-none_ME0.mat

%% SETUP METROPOLIS
options.mcmc.nite=100000;
options.mcmc.i_plot=20000;
n_reals_out=200;
options.mcmc.i_sample=options.mcmc.nite/n_reals_out;
randn('seed',2);rand('seed',2);


options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT PRIOR and POSTERIOR MOVIES
sippi_plot_movie(options.txt);

