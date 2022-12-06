clear all;close all
load caseBayesian_dx50_Ffat-none_ME0

N = 10000;
% sample prior
[m,prior]=sippi_prior(prior);
for i=1:N;
    if mod(i,100)==0, progress_txt(i,N,'Forward');end
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
T=100
Pacc = exp( (1/T)*(logL-max(logL)) );
r=rand(1,N);
i_sample = find(Pacc>r);

m_post=m_propose(:,:,i_sample);

[m_mean,m_var]  = etype(m_post);
figure;
subplot(1,2,1)
imagesc(prior{1}.x,prior{1}.y,m_mean);
caxis(cax);axis image
colorbar
subplot(1,2,2)
imagesc(prior{1}.x,prior{1}.y,sqrt(m_var));
colorbar
%caxis([0 0.005]);
axis image



save(sprintf('%s_rejection_out',txt))

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

