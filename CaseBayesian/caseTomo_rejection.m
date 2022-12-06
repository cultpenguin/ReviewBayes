clear all;close all
load caseBayesian_dx50_Ffat-none_ME0

nx=length(prior{1}.x);
ny=length(prior{1}.y);

N = 100000;
m_propose=zeros(ny,nx,N);
d_propose=zeros(length(d{1}),N);
logL=zeros(1,N);

% sample prior
[m,prior]=sippi_prior(prior);
t0=now;
parfor i=1:N;
    if mod(i,1000)==0, 
        [t_end_txt,t_left_seconds]=time_loop_end(t0,i,N);
        progress_txt(i,N,['Forward ',t_end_txt]);
    end
    % sample prior;
    m=sippi_prior(prior);
    m_propose(:,:,i)=m{1};
    % compute forward response
    d=sippi_forward(m,forward,prior);
    d_propose(:,i)=d{1};
    % compute log-likelihood
    logL(i)=sippi_likelihood(d,data);
end
%% simulate noise (for ML and EnK)%
doSimNoise = 1;
if doSimNoise==1
    Ct = diag(data{1}.d_std) + data{1}.Ct;
    try
        t0 = data{1}.t0;
    catch
        t0 = 0
    end
    d_noise=gaussian_simulation_cholesky(t0,Ct,N);
end
d_sim=d_noise.*0;
for i=1:N;
    d_sim(:,i)=d_propose(:,i)+d_noise(:,i);
end

%% Rejection
T=30
Pacc = exp( (1/T)*(logL-max(logL)) );
r=rand(1,N);
i_sample = find(Pacc>r);

m_post=m_propose(:,:,i_sample);
n_post=size(m_post,3)

[m_mean,m_var]  = etype(m_post);

%save(sprintf('%s_rejection_out',txt),'-v7.3')
%% SAVE SA HDF5??

%%

nx=length(prior{1}.x);
ny=length(prior{1}.y);

dx=prior{1}.x(2)-prior{1}.x(1);

Nr=5;

figure(11);
for i=1:Nr
    subplot(3,Nr,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,m_propose(:,:,i))
    axis image
    caxis(prior{1}.cax);colormap(cmap)
    title('\rho(m)\rightarrowm^*')
end
print_mul(sprintf('%s_N%d_prior_sample',txt,N))
figure(12);
for i=1:Nr
    subplot(3,Nr,i)
    try
    imagesc(prior{1}.x,prior{1}.y,m_post(:,:,i))
    end
    axis image;caxis(prior{1}.cax);colormap(cmap)
    title('\sigma(m)\rightarrowm^*')
end
print_mul(sprintf('%s_N%d_post_sample',txt,N))

figure(13);clf
subplot(1,2,1)
imagesc(prior{1}.x,prior{1}.y,m_mean)
axis image;caxis(prior{1}.cax);colormap(cmap)
colorbar
title('\sigma(m) - mean')
subplot(1,2,2)
imagesc(prior{1}.x,prior{1}.y,sqrt(m_var))
axis image;colormap(cmap)
title('\sigma(m) - standard deviation')
colorbar
print_mul(sprintf('%s_N%d_post_mean_std',txt,N))
    
save(sprintf('%s_out',txt))

