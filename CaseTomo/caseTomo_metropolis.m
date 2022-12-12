% clear all
close all

if ~exist('fmat','var'); fmat='caseBayesian_dx50_Ffat-none_ME0';end
if ~exist('N','var');
    N=10000;
    N = ceil(300000/32);
end
if ~exist('doSave','var'); doSave=1; end
if ~exist('rseed','var'); rseed=1;end
if rseed>0; rng('default') ;rng(rseed);end


%%
load(fmat,'prior','forward','data','txt')

%% SETUP METROPOLIS
options.mcmc.nite=N;
options.mcmc.i_plot=ceil(options.mcmc.nite/10);
options.mcmc.n_reals=200;
%i_sample=ceil(options.mcmc.nite/n_reals_out);
randn('seed',2);rand('seed',2);


options=sippi_metropolis(data,prior,forward,options);

%%
[reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(options.txt,1,100,1);

%%
if doSave==1
    save(sprintf('%s_metropolis_out',txt))
end

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
    caxis(prior{1}.cax);colormap(prior{1}.cmap)
    title('\rho(m)\rightarrowm^*')
end
print_mul(sprintf('%s_N%d_prior_sample',txt,options.mcmc.nite))
figure(12);
for i=1:Nr
    subplot(3,Nr,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,reals(:,:,i))
    axis image;caxis(prior{1}.cax);colormap(prior{1}.cmap)
    title('\sigma(m)\rightarrowm^*')
end
print_mul(sprintf('%s_N%d_post_sample',txt,options.mcmc.nite))

figure(13);clf
subplot(1,2,1)
imagesc(prior{1}.x,prior{1}.y,reshape(mean(reals_all),ny,nx))
axis image;caxis(prior{1}.cax);colormap(prior{1}.cmap)
colorbar
title('\sigma(m) - mean')
subplot(1,2,2)
imagesc(prior{1}.x,prior{1}.y,reshape(std(reals_all),ny,nx))
axis image;colormap(prior{1}.cmap)
title('\sigma(m) - standard deviation')
colorbar
print_mul(sprintf('%s_N%d_post_mean_std',txt,options.mcmc.nite))
    
%%
figure(15);clf
plot(options.C{1}.mcmc.logL)
xl=xlim;
xlim([xl(2)*.01 N])
ylim([-105 -50])
grid on
xlabel('Iteration Number')
ylabel('log(L(m))')
ppp(12,8,10,2,2)
print_mul(sprintf('%s_N%d_logL',txt,options.mcmc.nite))


