% caseTomo_lookup

% clear all
close all

if ~exist('fmat','var'); fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat';end
if ~exist('N','var');
    N=10000;
    N = ceil(300000/32);
end
if ~exist('Nme','var');Nme=0;end
if ~exist('doSave','var'); doSave=1; end
if ~exist('T_end','var'); T_end=1; end
if ~exist('rseed','var'); rseed=1;end
if rseed>0; rng('default') ;rng(rseed);end

%%
load(fmat,'prior','forward','data','txt')

txt_out = sprintf('%s_lu_N%d',txt,N);
disp(txt_out)

if exist([txt_out,'.mat'],'file')
    load(txt_out)
else
    
    txt_out_nn=sprintf('%s_lu',txt_out);
    if exist([txt_out_nn,'.mat'],'file')
        disp(sprintf('Loading %s',txt_out_nn))
        load(txt_out_nn,'ABC')
    else
        ABC=sippi_abc_setup(prior,forward,N,Nme);
        save(txt_out_nn,'ABC','-v7.3')
    end
    ns=200;
    [logL,evidence,T_est,ABC,dt,iCT]=sippi_abc_logl(ABC,data);
    [m_real, P_acc, i_use_all, d_real] = sippi_abc_post_sample(ABC, ns, T_est, logL);
    reals_all = m_real{1}';
end
%% sample prior
clear prior_reals
for i=1:10;
    m=sippi_prior(prior);
    prior_reals(i,:)=m{1}(:);
end

%% plot Post Stats
caseTomo_plot_post_stats(reals_all,prior_reals,prior,txt_out);
