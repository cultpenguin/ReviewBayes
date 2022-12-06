clear all;close all
%load caseBayesian_dx15_Fray_2d-none_ME0.mat
load caseBayesian_dx50_Fray_2d-none_ME0
%load caseBayesian_dx25_Feikonal-none_ME0.mat
%load caseBayesian_dx50_Feikonal-none_ME0.mat
load('caseBayesian_dx50_Ffat-none_ME0')

%% SETUP METROPOLIS
options.mcmc.nite=200000;
options.mcmc.i_plot=ceil(options.mcmc.nite/10);
n_reals_out=200;
options.mcmc.i_sample=options.mcmc.nite/n_reals_out;
randn('seed',2);rand('seed',2);


options=sippi_metropolis(data,prior,forward,options);

%%
[reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(options.txt,1,100,1);

%%
%% 

nx=length(prior{1}.x);
ny=length(prior{1}.y);

dx=prior{1}.x(2)-prior{1}.x(1);

Nr=5;

figure(11);
for i=1:Nr
    subplot(3,Nr,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,m{1})
    axis image
    caxis(prior{1}.cax);colormap(cmap)
    title('\rho(m)\rightarrowm^*')
end
print_mul(sprintf('%s_N%d_prior_sample',txt,options.mcmc.nite))
figure(12);
for i=1:Nr
    subplot(3,Nr,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,reals(:,:,i))
    axis image;caxis(prior{1}.cax);colormap(cmap)
    title('\sigma(m)\rightarrowm^*')
end
print_mul(sprintf('%s_N%d_post_sample',txt,options.mcmc.nite))

figure(13);clf
subplot(1,2,1)
imagesc(prior{1}.x,prior{1}.y,reshape(mean(reals_all),ny,nx))
axis image;caxis(prior{1}.cax);colormap(cmap)
colorbar
title('\sigma(m) - mean')
subplot(1,2,2)
imagesc(prior{1}.x,prior{1}.y,reshape(std(reals_all),ny,nx))
axis image;colormap(cmap)
title('\sigma(m) - standard deviation')
colorbar
print_mul(sprintf('%s_N%d_post_mean_std_N',txt,options.mcmc.nite))
    
save(sprintf('%s_out',txt))

return


%% PLOT SAMPLE FROM PRIOR
%sippi_plot_prior(options.txt);

%% PLOT SAMPLE FROM POSTERIOR
%sippi_plot_posterior(options.txt);

%% PLOT PRIOR and POSTERIOR MOVIES
%sippi_plot_movie(options.txt);

