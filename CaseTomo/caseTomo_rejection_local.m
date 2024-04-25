%clear all;
close all



if ~exist('fmat','var'); fmat='caseBayesian_dx50_Ffat-none_ME0';end
if ~exist('N','var');
    N=10000;
    N = ceil(300000/32);
end
if ~exist('doSave','var'); doSave=1; end
if ~exist('doSimNoise','var'); doSimNoise=1; end
if ~exist('rseed','var'); rseed=1;end
if ~exist('di_use','var'); di_use=1; end
if rseed>0; rng('default');rng(rseed);end


%%
load(fmat,'prior','forward','data','txt')
Nd=length(data{1}.d_obs);
if di_use>1
    data{1}.i_use = 1:di_use:length(data{1}.d_obs);
    Nd=length(data{1}.i_use);
end
nx=length(prior{1}.x);
ny=length(prior{1}.y);

txt_h5=sprintf('%s_rejection_N%d_di%d_out',txt,N,di_use)
h5=[txt_h5,'.h5'];

if exist([txt_h5,'.mat'],'file')
    load(txt_h5)
else
    %%
    clc
    delete(h5)
    h5create(h5,'/D',[Nd,N],'ChunkSize',[Nd,1])
    h5create(h5,'/Dsim',[Nd,N],'ChunkSize',[Nd,1])
    h5create(h5,'/M',[ny,nx,N],'ChunkSize',[ny nx 1])

    if doSimNoise==1
        try
            Ct = data{1}.Ct;
        catch
            Ct = diag(data{1}.d_std);
        end
        try
            dt = data{1}.dt;
        catch
            dt = 0;
        end
        % get choleskey
        [~,~,Ct_chol]=gaussian_simulation_cholesky(dt,Ct,1);
    end

    %%
    % initialize
    [m,prior]=sippi_prior(prior);
    [d,forward]=sippi_forward(m,forward,prior,data);
    [~,~,data]=sippi_likelihood(d,data);

    N_chunk = 250000;
    N_chunk = 100000;
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
        d_sim=rand(Nd,N_chunk);
        


        Nd_in_loop=(1+i_end-i_start);
        %for i=1:Nd_in_loop
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
            d=sippi_forward(m,forward,prior,data);
            d_propose(:,i)=d{1};
            % compute log-likelihood
            logL_propose(i)=sippi_likelihood(d,data);

            if doSimNoise==1
                %d_noise=gaussian_simulation_cholesky(dt,Ct,1);
                is_cholesky=1;
                d_noise=gaussian_simulation_cholesky(dt,Ct_chol,1,is_cholesky);
                if di_use>1
                    d_noise=d_noise(data{1}.i_use);
                end
                d_sim(:,i)=d_propose(:,i)+d_noise(:);
            end            
        end
        t_end = now;
        t_elapsed_minutes=(t_end-t_start)*60*24;
        disp(sprintf('Time to setup [M,D]_[%d/%d]: %4.3f minutes',i_end,N,t_elapsed_minutes))
        logL(i_start:1:i_end)=logL_propose(1:Nd_in_loop);
        h5write(h5,'/M',m_propose(:,:,1:Nd_in_loop),[1,1,i_start],[ny,nx,Nd_in_loop])
        h5write(h5,'/D',d_propose(:,1:Nd_in_loop),[1,i_start],[Nd,Nd_in_loop])
        h5write(h5,'/Dsim',d_sim(:,1:Nd_in_loop),[1,i_start],[Nd,Nd_in_loop])

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

Pacc = exp( (1./T)*(logL-max(logL)) );

r=rand(1,N);
i_sample = find(Pacc>r);

%m_propose_big = h5read(h5,'/M');
%m_post_big=m_propose_big(:,:,i_sample);
for i=1:length(i_sample)
    m_post(:,:,i)=h5read(h5,'/M',[1 1 i_sample(i)],[ny nx 1]);
    mm=m_post(:,:,i);
    post_reals(i,:)=mm(:);
    post_reals_2d(:,:,i)=mm;
end
n_post=size(m_post,3);

[m_mean,m_var]  = etype(m_post);

txt_out = sprintf('%s_rejection_N%d_di%d_T%d',txt,N,di_use,T);
disp(txt_out)

post_mean = reshape(m_mean,ny,nx);
post_std = sqrt(reshape(m_var,ny,nx));

h5writeMatrix(h5,'/post_reals',post_reals_2d)
h5writeMatrix(h5,'/post_mean',post_mean)
h5writeMatrix(h5,'/post_std',post_std)
h5writeMatrix(h5,'/T',T)


%& save post_reals, post_mean, post_std

%%

%% sample prior
clear prior_reals
for i=1:10;
    m=sippi_prior(prior);
    prior_reals(i,:)=m{1}(:);
end

%% Convert to velocity?
if forward.is_slowness == 1
    prior_reals = 1./prior_reals;
    post_reals = 1./post_reals;
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



%% LOCAL
close all
M=h5read(h5,'/M');
D=h5read(h5,'/D');
Dsim=h5read(h5,'/Dsim');
forward_lin=load('caseTomo_Kallerup_dx10_Ffat-none_ME0_slo1_SE2_G0.mat','forward');
G=forward_lin.forward.G;
%%
dx = prior{1}.x(2)-prior{1}.x(1);
local_x = 1;local_y = local_x;
local_x = 1.25;local_y = 1.75;
%local_x = .5;local_y = .7;
%local_x = 2.4;local_y = local_x;
wx = ceil(local_x/dx);
wy = ceil(local_y/dx);
%wy = wx;
x1 = 0:wx:length(prior{1}.x);
y1 = 0:wy:length(prior{1}.y);
if length(prior{1}.x)>(x1(end)+1);
    x1 = [x1,length(prior{1}.x)];
end
if length(prior{1}.y)>(y1(end)+1);
    y1 = [y1,length(prior{1}.y)];
end
m=sippi_prior(prior);
m0 = m{1}.*0;

ns=400;
post_mean_local = m0;
post_std_local = m0;
post_reals_local=zeros(ny,nx,ns);
clear T_local
clear t
clear ndata 
t0=now;
k=0;
for ix = 1:(length(x1)-1)
for iy = 1:(length(y1)-1)
    k=k+1;
    x = int16((x1(ix)+1):1:x1(ix+1));
    y = int16((y1(iy)+1):1:y1(iy+1));
    [xx,yy]=meshgrid(x,y);
    m_use = m{1}.*0;
    m_use(yy(:),xx(:))=1;
    ixy_use = find(m_use(:)>0);
    dg=zeros(1,size(G,1));
    for ig =1:size(G,1);
         dg(ig)=sum(abs(G(ig,ixy_use)));
    end
    id_use= find(dg>=0);
    id_use= find(dg>0.001);
    %id_use= find(dg>0.1);
    %id_use= find(dg>0.05);
    id_use= find(dg>0.5);
    ndata(k)=length(id_use);
    iCD=inv(data{1}.CD(id_use,id_use));
    logL=zeros(1,N);
    for i=1:N
        dd=(data{1}.d_obs-D(:,i))-data{1}.dt;
        dd=dd(id_use);
        logL(i)=-.5*dd'*iCD*dd;
    end
    sort_logL=sort(logL-max(logL));
    if sum(~isnan(sort_logL))>0
        if isnan(sort_logL(end))
            sort_logL=sort_logL(~isnan(sort_logL));
        end
        logl_lev=sort_logL(end-N_above-1);
        
        T_est_local = logl_lev/log(P_acc_lev);
        T_est_local = max([1 T_est_local]);
    else
        T_est_local = inf;
    end
    T_est_local
    %
    
    T_local(k)=1;
    T_local(k) = T_est_local
    [i_use_all,P_acc]=sippi_abc_post_sample_logl(logL,ns,T_local(k));

    m_post = M(:,:,i_use_all);
    m_post_mean = mean(m_post,3);
    m_post_std = m_post_mean*0;
    
    
    for jx=1:nx
    for jy=1:ny
        m_post_std(jy,jx)=std(squeeze(m_post(jy,jx,:)));
    end
    end
    
    for jx=1:nx
    for jy=1:ny
        if m_use(jy,jx)==1
            post_reals_local(jy,jx,:)=m_post(jy,jx,:);
            post_mean_local(jy,jx)=m_post_mean(jy,jx);
            post_std_local(jy,jx)=m_post_std(jy,jx);
        end
    end
    end

    t(k)=(now-t0)*3600*24;

    subplot(1,4,1)
    imagesc(prior{1}.x,prior{1}.y,reshape(sum(G(id_use,:),1),ny,nx))
    axis image;
    subplot(1,4,2)
    imagesc(prior{1}.x,prior{1}.y,m_use)
    axis image
    subplot(1,4,3)
    imagesc(prior{1}.x,prior{1}.y,post_mean_local)
    caxis([0.0700    0.1800])
    colormap(gca,jet)
    axis image
    subplot(1,4,4)
    imagesc(prior{1}.x,prior{1}.y,post_std_local)
    axis image
    
    drawnow
    pause(.01)
    
end
end

for i=1:ns
    mm=post_reals_local(:,:,i);
    post_reals_local_flat(i,:)=mm(:);
end


%& save post_reals, post_mean, post_std
txt_h5_local = sprintf('%s_localized_wx%d_wy%d_T%d',txt_h5,wx,wy,T);
copyfile([txt_h5,'.h5'],[txt_h5_local,'.h5'])
h5writeMatrix(h5,'/post_reals_local',post_reals_local)
h5writeMatrix(h5,'/post_mean_local',post_mean_local)
h5writeMatrix(h5,'/post_std_local',post_std_local)
try;h5writeMatrix(h5,'/T_local',T_local);end
h5writeMatrix(h5,'/window',[wx,wy])


if doSave==1
    save([txt_h5_local,'.mat'],'post_*','T*','wx','wy','t','local_x','local_y','ndata')
end

%%
caseTomo_plot_post_stats(post_reals_local_flat,prior_reals,prior,[txt_out,'_local']);
