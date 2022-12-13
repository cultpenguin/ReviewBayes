%clear all;
close all

if ~exist('fmat','var'); fmat='caseBayesian_dx50_Ffat-none_ME0';end
if ~exist('N','var');
    N=10000;
    N = ceil(300000/32);
end
if ~exist('doSave','var'); doSave=1; end
if ~exist('doSimNoise','var'); doSimNoise=0; end
if ~exist('rseed','var'); rseed=1;end
if rseed>0; rng('default');rng(rseed);end


%%
load(fmat,'prior','forward','data','txt')

nx=length(prior{1}.x);
ny=length(prior{1}.y);
Nd=length(data{1}.d_obs);

txt_h5=sprintf('%s_rejection_N%d_out',txt,N)
h5=[txt_h5,'.h5'];

if exist([txt_h5,'.mat'],'file')
    load(txt_h5)
else
    %%
    clc
    delete(h5)
    h5create(h5,'/D',[Nd,N],'ChunkSize',[Nd,1])
    h5create(h5,'/M',[ny,nx,N],'ChunkSize',[ny nx 1])


    %%
    % initialize
    [m,prior]=sippi_prior(prior);
    [d,forward]=sippi_forward(m,forward,prior);
    [~,~,data]=sippi_likelihood(d,data);

    N_chunk = 50000;
    N_loop = ceil(N/N_chunk);

    logL=zeros(1,N);
    t_start=now;
    for ic = 1:N_loop
        i_start = (ic-1)*N_chunk+1;
        if ic==N_loop
            i_end = N;
        else
            i_end = ic*N_chunk;
        end
        disp(sprintf('ic=%02d/%02d, %06d-%06d -- N_chunk=%d',ic,N_loop,i_start,i_end,N_chunk))
        ii=i_start:1:i_end;
        m_propose=rand(ny,nx,N_chunk);
        d_propose=rand(Nd,N_chunk);
        logL_propose=zeros(1,N_chunk);



        Nd_in_loop=(1+i_end-i_start);
        parfor i=1:Nd_in_loop
            if mod(i,100)==0,
                [t_end_txt,t_left_seconds]=time_loop_end(t_start,i,Nd_in_loop);
                [t_end_txt,t_left_seconds]=time_loop_end(t_start,ii(i),N);
                progress_txt(ii(i),N,['Forward ',t_end_txt]);
            end
            % sample prior;
            m=sippi_prior(prior);
            m_propose(:,:,i)=m{1};
            % compute forward response
            d=sippi_forward(m,forward,prior);
            d_propose(:,i)=d{1};
            % compute log-likelihood
            logL_propose(i)=sippi_likelihood(d,data);
        end
        t_end = now;
        t_elapsed_minutes=(t_end-t_start)*60*24;

        logL(i_start:1:i_end)=logL_propose(1:Nd_in_loop);
        h5write(h5,'/M',m_propose(:,:,1:Nd_in_loop),[1,1,i_start],[ny,nx,Nd_in_loop])
        h5write(h5,'/D',d_propose(:,1:Nd_in_loop),[1,i_start],[Nd,Nd_in_loop])

    end

    clear m_propose
    clear d_propose


    %logL(i_start:1:i_end)=logL_propose(1:Nd_in_loop);
    %h5write(h5,'/M',m_propose(:,:,1:Nd_in_loop),[1,1,i_start],[ny,nx,Nd_in_loop])
    %h5write(h5,'/D',d_propose(:,1:Nd_in_loop),[1,i_start],[Nd,Nd_in_loop])

    disp(sprintf('Time to setup [M,D]: %4.3f minutes',t_elapsed_minutes))


    %%
    %d5=h5read(h5,'/D');
    %m5=h5read(h5,'/M');
    %% simulate noise (for ML and EnK)%
    if doSimNoise==1
        Ct = diag(data{1}.d_std) + data{1}.Ct;
        try
            t0 = data{1}.t0;
        catch
            t0 = 0
        end
        d_noise=gaussian_simulation_cholesky(t0,Ct,N);
        d_sim=d_noise.*0;
        for i=1:N;
            d_sim(:,i)=d_propose(:,i)+d_noise(:,i);
        end
    end

    if doSave==1
        save(txt_h5)
    end

end
%% SAVE AS HDF5??

%% Rejection
% Compute evidence
maxlogL=max(logL);
log_evidence = maxlogL + log( nansum(exp(logL-maxlogL))/length(logL) );
disp(sprintf('%s: logEvidence=%5.3f',txt_h5,log_evidence))

% Compute annealing temperature
N_above=20;
P_acc_lev=0.1;

sort_logL=sort(logL-max(logL));
if sum(~isnan(sort_logL))>0
    if isnan(sort_logL(end))
        sort_logL=sort_logL(~isnan(sort_logL));
    end
    logl_lev=sort_logL(end-N_above-1);

    T_est = logl_lev/log(P_acc_lev);
    T_est = max([1 T_est]);
else
    T_est = inf;
end
%
T=ceil(T_est);

Pacc = exp( (1/T)*(logL-max(logL)) );

r=rand(1,N);
i_sample = find(Pacc>r);

%m_propose_big = h5read(h5,'/M');
%m_post_big=m_propose_big(:,:,i_sample);
for i=1:length(i_sample)
    m_post(:,:,i)=h5read(h5,'/M',[1 1 i_sample(i)],[ny nx 1]);
    mm=m_post(:,:,i);
    post_reals(i,:)=mm(:);
end
n_post=size(m_post,3);

[m_mean,m_var]  = etype(m_post);

txt_out = sprintf('%s_rejection_N%d_T%d',txt,N,T);
disp(txt_out)

%%

%% sample prior
clear prior_reals
for i=1:10;
    m=sippi_prior(prior);
    prior_reals(i,:)=m{1}(:);
end

%% plot Post Stats
caseTomo_plot_post_stats(post_reals,prior_reals,prior,txt_out);
%%
figure(10);clf;
subplot(3,1,1)
semilogy(exp((logL-max(logL))),'.');ylim([1e-100 1])
title(sprintf('T=%g',1))
subplot(3,1,2)
semilogy(exp((1/T)*(logL-max(logL))),'.');ylim([1e-100 1])
title(sprintf('T=%g',T))
subplot(3,1,3)
semilogy(exp((1/T)*(logL-max(logL))),'-');ylim([1e-2 1])
title(sprintf('T=%g',T))
print_mul(sprintf('%s_anneal',txt_out))

return

%%
nx=length(prior{1}.x);
ny=length(prior{1}.y);

dx=prior{1}.x(2)-prior{1}.x(1);

Nr=25;
Nr_show=5;

figure(11);
for i=1:Nr
    %m=sippi_prior(prior);
    m{1}=h5read(h5,'/M',[1 1 i],[ny nx 1]);
    d=sippi_forward(m,forward,prior);
    d_prior(:,i)=d{1};
    if i<=Nr_show
        subplot(3,Nr_show,i)
        imagesc(prior{1}.x,prior{1}.y,m{1})
        axis image
        caxis(prior{1}.cax);colormap(prior{1}.cmap)
        title('\rho(m)\rightarrowm^*')
    end
end
print_mul(sprintf('%s_prior_sample',txt_out))
figure(12);
for i=1:Nr
    try
        m{1}=m_post(:,:,i);
        d=sippi_forward(m,forward,prior);
        d_post(:,i)=d{1};
        if i<=Nr_show
            subplot(3,Nr_show,i)
            imagesc(prior{1}.x,prior{1}.y,m{1})
            axis image;caxis(prior{1}.cax);colormap(prior{1}.cmap)
            title('\sigma(m)\rightarrowm^*')
        end
    end

end
print_mul(sprintf('%s_post_sample',txt_out))

figure(13);clf
subplot(1,2,1)
imagesc(prior{1}.x,prior{1}.y,m_mean)
axis image;caxis(prior{1}.cax);colormap(prior{1}.cmap)
colorbar
title('\sigma(m) - mean')
subplot(1,2,2)
imagesc(prior{1}.x,prior{1}.y,sqrt(m_var))
axis image;colormap(prior{1}.cmap)
%colormap(gca,cmap_linear([1 1 1;0 0 0]))
try
    caxis(prior{1}.cax_std);
catch
    caxis([0 0.02]);
end
title('\sigma(m) - standard deviation')
colorbar
print_mul(sprintf('%s_post_mean_std',txt_out))

%%
figure(14);set_paper('landscape')
%try
%    d_std = sqrt(diag(data{1}.CD));
%catch
d_std = data{1}.d_std;
%end
p1=plot(d_prior,'k-','LineWidth',4,'color',[1 1 1].*.5);
hold on
pe=errorbar(data{1}.d_obs,2*d_std,'color','k');
p2=plot(d_post,'r-','LineWidth',.1);
hold off
grid on
xlim([1 702])
ylim([25 60])
print_mul(sprintf('%s_datafit',txt_out))




