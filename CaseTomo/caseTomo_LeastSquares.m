clear all;close all


% load data
%load('caseBayesian_dx50_Fray_2d-none_ME0.mat','prior','data','D','cax','txt')
load('caseBayesian_dx10_Ffat-none_ME0','prior','data','D','cax','txt')
%load('caseBayesian_dx25_Ffat-none_ME0','prior','data','D','cax','txt')
%load('caseBayesian_dx50_Ffat-none_ME0','prior','data','D','cax','txt')


if ~exist('cmap','var'); cmap=jet;;end
if ~exist('rseed','var'); rseed=1;end
if rseed>0; rng('default') ;rng(rseed);end

%% FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.forward_function='sippi_forward_traveltime';
forward.is_slowness=1; % use slowness parameterization

%% THE PRIOR

if forward.is_slowness==1
    % slowness
    prior_org=prior;
    clear prior
    d = 1./prior_org{1}.o_nscore.d;
    m0 = mean(d);
    var0 = var(d);
    Va=deformat_variogram(prior_org{1}.Va);
    Va.par1=var0;
    prior{1}.Va=Va;
    prior{1}.m0=m0;
    prior{1}.name=prior_org{1}.name;
    prior{1}.type=prior_org{1}.type;
    prior{1}.x=prior_org{1}.x;
    prior{1}.y=prior_org{1}.y;
    prior{1}.cax=prior_org{1}.cax;
    %prior{1}=rmfield(prior{1},'fftma_options');
    %prior{1}=rmfield(prior{1},'seq_gibbs');
    %prior{1}=rmfield(prior{1},'o_nscore');
    %[d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
    %prior=sippi_prior_init(prior)
    %prior{1}.o_nscore=o_nscore;
end

nx=length(prior{1}.x);
ny=length(prior{1}.y);


m=sippi_prior(prior);

%% SOLVE LEAST SQUARES OPTIONS
options.lsq.type='lsq';
%lsq_type='visim';
%lsq_type='error_sim';

% set number of realization
options.lsq.n_reals=150;

% select a subset of data to consider
%data{1}.i_use=1:20:702;



%% STRAIGHT RAY FORWARD
forward.forward_function='sippi_forward_traveltime';
forward.type='ray';
forward.linear=1;
options.txt=[txt,'_',forward.type];
[m_est_1,Cm_est_1,m_reals_1,options_1,data_1,prior_1,forward_1]=sippi_least_squares(data,prior,forward,options);

%% LINEAR FAT FORWARD
try; forward=rmfield(forward,'G');end
forward.type='fat';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[m_est_2,Cm_est_2,m_reals_2,options_2,data_2,prior_2,forward_2]=sippi_least_squares(data,prior,forward,options);



%% 
dx=prior{1}.x(2)-prior{1}.x(1);
txt=sprintf('caseTomo_lsq_%s_dx%d',forward.type,ceil(100*dx));
Nr=5;

figure(11);
for i=1:Nr
    subplot(3,Nr,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m{1},ny,nx))
    axis image
    caxis(prior{1}.cax);colormap(cmap)
    title('\rho(m)\rightarrowm^*')
end
print_mul(sprintf('%s_prior_sample',txt))

figure(12);
for i=1:Nr
    subplot(3,Nr,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_1(:,i),ny,nx))
    %imagesc(prior{1}.x,prior{1}.y,reshape(1./m{1},ny,nx))
    axis image;caxis(prior{1}.cax);colormap(cmap)
    title('\sigma(m)\rightarrowm^*')
end
print_mul(sprintf('%s_post_sample',txt))

figure(13);clf
subplot(1,2,1)
imagesc(prior{1}.x,prior{1}.y,reshape(mean(1./m_reals_1'),ny,nx))
axis image;caxis(prior{1}.cax);colormap(cmap)
colorbar
title('\sigma(m) - mean')
subplot(1,2,2)
imagesc(prior{1}.x,prior{1}.y,reshape(std(1./m_reals_1'),ny,nx))
axis image;colormap(cmap)
title('\sigma(m) - standard deviation')
colorbar
print_mul(sprintf('%s_post_mean_std',txt))
    


return


%% POST PLOTS
%sippi_plot_posterior_sample(options_1.txt)
%sippi_plot_posterior_sample(options_2.txt)
%sippi_plot_posterior_sample(options_3.txt)

%%
cax=fliplr(1./prior{1}.cax);
figure(3);clf;
if useSynth==1;
    subplot(1,4,1);
    imagesc(prior{1}.x,prior{1}.y,1./m_ref{1});caxis(cax);axis image
    title('Reference')
end
subplot(1,4,2);
imagesc(prior{1}.x,prior{1}.y,1./m_est_1{1});caxis(cax);axis image
title('a) Ray kernel')
subplot(1,4,3);
imagesc(prior{1}.x,prior{1}.y,1./m_est_2{1});caxis(cax);axis image
title('a) Fat kernel')
subplot(1,4,4);
imagesc(prior{1}.x,prior{1}.y,1./m_est_3{1});caxis(cax);axis image
colorbar_shift;
title('a) Born kernel')
print_mul(sprintf('%s_compare_est',txt));

%%
figure(4);clf,set_paper('landscape')
for i=1:5;
    
    subplot(3,5,i);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_1(:,i),ny,nx))
    axis image;caxis(cax);
    
    %subplot(3,5,i+5);
    %imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_2(:,i),ny,nx))
    %axis image;caxis(cax);
    
    %subplot(3,5,i+10);
    %imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_3(:,i),ny,nx))
    %axis image;caxis(cax);
    
end
subplot(3,5,10);colorbar_shift;
try
    n_use=length(data{1}.i_use);
catch
    n_use=length(data{1}.d_obs);
end
print_mul(sprintf('%s_compare_reals_nd%d',txt,n_use));

%% 
Nr=5;
for i=1:Nr
    subplot(3,5,i)
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m{1},ny,nx))
    axis image
    title('\rho(m)\rightarrowm^*')
end
