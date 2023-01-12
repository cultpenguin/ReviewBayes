% clear all
close all

if ~exist('fmat','var'); fmat='caseBayesian_dx50_Ffat-none_ME0';end
if ~exist('N','var');
    N=10000;
    N = ceil(300000/32);
end
if ~exist('doSave','var'); doSave=1; end
if ~exist('T_end','var'); T_end=1; end
if ~exist('rseed','var'); rseed=1;end
if rseed>0; rng('default') ;rng(rseed);end

%%
load(fmat,'prior','forward','data','txt')

txt_out = sprintf('%s_metropolis_N%d',txt,N);
disp(txt_out)

if exist([txt_out,'.mat'],'file')
    load(txt_out)
else

    %% SETUP METROPOLIS
    options.mcmc.nite=N;
    n1=ceil(N/10);
    options.mcmc.i_plot=ceil(options.mcmc.nite/10);
    options.mcmc.i_plot=ceil(options.mcmc.nite/50);
    options.mcmc.n_reals=200;
    %i_sample=ceil(options.mcmc.nite/n_reals_out);
    randn('seed',2);rand('seed',2);
    %options.mcmc.T=1;
    for ip=1:length(prior)
        prior{ip}.seq_gibbs.i_update_step_max=2*n1;
    end
    % ANNEALING (TEMPERATURE AS A FUNCTION OF ITERATION NUMBER)
    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=n1; %  iteration number when annealing stops
    options.mcmc.anneal.T_begin=25; % Start temperature for annealing
    options.mcmc.anneal.T_end=T_end; % End temperature for annealing
    
    try;options.txt=txt;end;
    options=sippi_metropolis(data,prior,forward,options);

    %%
    [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(options.txt,1,100,1);

    %%
    if doSave==1
        save(txt_out)
    end

end

%% sample prior
clear prior_reals
for i=1:10;
    m=sippi_prior(prior);
    prior_reals(i,:)=m{1}(:);
end

%% plot Post Stats
caseTomo_plot_post_stats(reals_all,prior_reals,prior,txt_out);

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
print_mul(sprintf('%s_logL',txt_out))


