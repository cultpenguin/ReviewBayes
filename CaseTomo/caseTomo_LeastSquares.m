close all

if ~exist('fmat','var'); fmat='caseBayesian_dx50_Ffat-none_ME0';end
if ~exist('rseed','var'); rseed=1;end
if ~exist('di_use','var'); di_use=1; end
if rseed>0; rng('default') ;rng(rseed);end

% load data
load(fmat,'prior','data','forward','D','cax','txt')
if di_use>1
    data{1}.i_use = 1:di_use:length(data{1}.d_obs);
end

%if ~exist('cmap','var'); cmap=jet;;end

%% FORWARD MODEL
%forward.sources=D.S;
%forward.receivers=D.R;
%forward.forward_function='sippi_forward_traveltime';
%forward.is_slowness=1; % use slowness parameterization

%% Make sample of initial prior
clear prior_reals
for i=1:25;
    m=sippi_prior(prior);
    prior_reals_org(i,:)=m{1}(:);
end
if forward.is_slowness==0
    prior_reals_org=1./prior_reals_org;
end

%% THE PRIOR
d = prior{1}.o_nscore.d; % CONVERT TO GAUSSIAN PRIOR!!!
if forward.is_slowness==0
    d = 1./d;
    prior{1}.cax=fliplr(1./prior{1}.cax);
    forward.is_slowness=1;
    forward.linear_m=1./forward.linear_m;
end
prior{1}=rmfield(prior{1},'d_target');
prior{1}=rmfield(prior{1},'o_nscore');
m0 = mean(d);
var0 = var(d);
Va=deformat_variogram(prior{1}.Va);
Va.par1=var0;
prior{1}.Va=Va;
prior{1}.m0=m0;

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



%% LEAST SQUARES INVERSION
%forward.forward_function='sippi_forward_traveltime';
%forward.type='ray';
%forward.linear=1;
options.txt=[txt,'_',forward.type];
options.lsq.n_reals=500;
[m_est,Cm_est,post_reals,options_1,data_1,prior_1,forward_1]=sippi_least_squares(data,prior,forward,options);

clear prior_reals
for i=1:25;
    m=sippi_prior(prior);
    prior_reals(i,:)=m{1}(:);
end

%%
txt_out=sprintf('%s_di%d_LSQ',txt,di_use);

caseTomo_plot_post_stats(1./post_reals',1./prior_reals,prior,txt_out)


%%
figure
subplot(1,2,1)
histogram(prior_reals_org(:))
hold on
histogram(prior_reals(:))
hold off

subplot(1,2,2)
histogram(1./prior_reals_org(:))
hold on
histogram(1./prior_reals(:))
hold off
print_mul(sprintf('%s_hist_compare',txt_out))




