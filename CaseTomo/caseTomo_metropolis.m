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

txt_out = sprintf('%s_metropolis_N%d',txt,N);
disp(txt_out)

if exist([txt_out,'.mat'],'file')
    load(txt_out)
else

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


