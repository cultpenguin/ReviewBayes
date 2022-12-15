close all

if ~exist('fmat','var'); fmat='caseBayesian_dx50_Ffat-none_ME0';end
if ~exist('rseed','var'); rseed=1;end
if rseed>0; rng('default') ;rng(rseed);end

% load data
load(fmat,'prior','data','forward','D','cax','txt')

%if ~exist('cmap','var'); cmap=jet;;end

%% FORWARD MODEL
%forward.sources=D.S;
%forward.receivers=D.R;
%forward.forward_function='sippi_forward_traveltime';
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
    prior{1}.cmap=prior_org{1}.cmap;
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



%% LEAST SQUARES INVERSION
%forward.forward_function='sippi_forward_traveltime';
%forward.type='ray';
%forward.linear=1;
options.txt=[txt,'_',forward.type];
options.lsq.n_reals=500;
[m_est,Cm_est,post_reals,options_1,data_1,prior_1,forward_1]=sippi_least_squares(data,prior,forward,options);

clear prior_reals
for i=1:10;
    m=sippi_prior(prior);
    prior_reals(i,:)=m{1}(:);
end

%%
txt_out=sprintf('%s_LSQ',txt);
caseTomo_plot_post_stats(1./post_reals',1./prior_reals,prior_org,txt_out)

