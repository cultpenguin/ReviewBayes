%%
close all;
doPlot=1;

if ~exist('Nlu','var')
    Nlu=1000000;
    Nlu=500000;
    Nlu=100000;
    %Nlu=10000;
    %Nlu=1000;
end
if ~exist('Nr','var')
    Nr=400;
end
if ~exist('Cd','var')
    %Cd = [0.0025, -0.005; -.005, 0.0225]; % (r0,g)
    Cd = [0.0025, -0.005 0 ; -.005, 0.0225 0 ; 0 0 1^2]; % (r0,g)
    % Lower the noise, as the noise is probably too high!!!
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.125*sqrt(abs(Cd(1:2,1:2)))).^2;
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.25*sqrt(abs(Cd(1:2,1:2)))).^2;
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.5*sqrt(abs(Cd(1:2,1:2)))).^2;
end

%% SIMPLE TEST
doSimpleTest=0;
if doSimpleTest==1;
    D1=load('data_cube_UPDATED.mat');  Nlu=Nlu+1;
    D2=load('data_cube_UPDATED2.mat');    Nlu=Nlu+2;
    D2=load('data_cube_UPDATED3_no_corr');   Nlu=Nlu+3;
    d=repmat(D1.depth,[100,1]);
    ii=find(d>2140 & d<2141 );
    %ii=find(d>2100 & d<2101 );
    %ii=find(d>2200 & d<2201 );
    %ii=25:25:4414400;
    %ii=1:10:4414400/100;
    subplot(1,2,1);
    scatter(D1.g_sim_no_noise(ii), D1.r0_sim_no_noise(ii), 15, D1.v_p(ii), 'filled')
    colorbar;%caxis([2000 3000])
    %hold on;plot(D1.g_data,D1.r_0_data,'k.');hold off
    axis([-.5 .5 -.3 .3])
    ax=axis;
    subplot(1,2,2);
    scatter(D2.g_sim_no_noise(ii), D2.r0_sim_no_noise(ii), 5, D2.v_p(ii), 'filled');
    colorbar;%caxis([2000 3000])
    %hold on;plot(D2.g_data,D2.r_0_data,'k.');hold off
    axis(ax)
    colormap jet
    return
end
%%
disp('Load data')
%D=load('data_cube_ALL.mat');
%D=load('data_cube_UPDATED.mat');  Nlu=Nlu+1;
%D=load('data_cube_UPDATED2.mat');    Nlu=Nlu+2;
D=load('data_cube_UPDATED3_no_corr');   Nlu=Nlu+3;

D.depth_all=repmat(D.depth,[100,1]);
ii=find(D.depth_all>2140 & D.depth_all<2141 );
x=unique(D.xline_obs);nx=length(x);
y=unique(D.inline_obs);ny=length(y);


%%
%
N = length(D.v_p);
ii = ceil(linspace(1,N,20000));

figure(1)
subplot(2,2,[1,2])
scatter(D.xline_obs,D.inline_obs,5,D.depth,'filled')
axis image;colorbar
xlabel('Xline');ylabel('Inline');title('depth')
subplot(2,2,3)
scatter(D.xline_obs,D.inline_obs,5,D.g_data,'filled')
axis image;colorbar
xlabel('Xline');ylabel('Inline');title('g')
subplot(2,2,4)
scatter(D.xline_obs,D.inline_obs,5,D.r_0_data,'filled')
axis image;colorbar
xlabel('Xline');ylabel('Inline');title('R0')
print_mul('Data')


figure(2);
subplot(2,1,1)
scatter(D.v_p(ii),D.rho(ii),4,D.g_sim(ii),'filled')
xlabel('V_p');ylabel('rho');title('g')
caxis([-1 1].*.3);colormap(gca,cmap_linear)
colorbar
subplot(2,1,2)
scatter(D.v_p(ii),D.rho(ii),4,D.r_0_sim(ii),'filled')
xlabel('V_p');ylabel('rho');title('R_0')
caxis([-1 1].*.2);colormap(gca,cmap_linear)
colorbar

%% Setup ABC
progress_out(sprintf('%s: Setup ABC start',mfilename))
disp('Setup ABC')
clear ABC
Nlu=min([N Nlu]);
i_use=randomsample(1:N,Nlu);
mat_file = sprintf('MINA_ABC_N%d',Nlu);

if exist([mat_file,'.mat'])==2
    disp(sprintf('Loading %s',mat_file))
    load(mat_file)
else

    i=0;
    for k=i_use;1:1000000;N;5000;%N;
        i=i+1;
        if mod(i,1000)==0;progress_txt(i,length(i_use));end
        % m
        useMstyle=1;
        if useMstyle==1;
            j=0;
            j=j+1;ABC.m{i}{j}=D.sat_g_log(k);
            j=j+1;ABC.m{i}{j}=D.sat_o_log(k);
            j=j+1;ABC.m{i}{j}=D.v_clay_log(k);
            j=j+1;ABC.m{i}{j}=D.v_p(k);
            j=j+1;ABC.m{i}{j}=D.v_s(k);
            j=j+1;ABC.m{i}{j}=D.rho(k);
            j=j+1;ABC.m{i}{j}=D.depth_all(k);
            %j=j+1;ABC.m{i}{j}=D.porosity(k);
        else
            ABC.m{i}{1}=[D.sat_g_log(k) D.sat_o_log(k) D.v_clay_log(k) D.v_p(k) D.v_s(k) D.rho(k) D.depth_all(k)];
        end

        if i==1;
            ABC.m=repmat(ABC.m,[length(i_use),1]);
        end

        % d
        j=0;
        %j=j+1;ABC.d{i}{j}=[D.g_sim(k), D.r_0_sim(k)];
        %j=j+1;ABC.d{i}{j}=[D.depth_all(k)];
        %j=j+1;ABC.d{i}{j}=[ D.r_0_sim(k), D.g_sim(k), D.depth_all(k)]';
        j=j+1;ABC.d{i}{j}=[ D.r0_sim_no_noise(k), D.g_sim_no_noise(k), D.depth_all(k)]';

        if i==1;
            ABC.d=repmat(ABC.d,[length(i_use),1]);
        end

    end
    disp(sprintf('Saving %s',mat_file))
    save(mat_file,'ABC');
end

Nlu=length(ABC.m);

progress_out(sprintf('%s: Setup ABC end',mfilename))
%%


%% test single inversion
ns=Nr;
clear d
id=1;
data=data_mul{id};

d{1}=data{1}.d_obs + gaussian_simulation_cholesky(0,data{1}.Cd);

[logL]=sippi_likelihood(d,data);

%% TEST SINGLE DATA
%ABC.use_sippi_likelihood=0;
id=1000;
data=data_mul{id};

t0=now;
[logL,evidence,T_est,ABC,dt,iCT]=sippi_abc_logl(ABC,data);
T=T_est
[m_post, P_acc, i_use_all,d_post] = sippi_abc_post_sample(ABC, ns, T, logL);

data_p=data;
data_p{1}.Cd=diag([100000,100000,1]);
[logL_p,evidence_p,T_p]=sippi_abc_logl(ABC,data_p);
[m_prior, P_acc_prior, i_use_all_prior,d_prior] = sippi_abc_post_sample(ABC, ns,T_p, logL_p);

t1=now;
t_per_data = (t1-t0)*3600*24;
t_total = t_per_data*Nd;
disp(sprintf('TimePerData=%4.1f ms',1000*t_per_data))
disp(sprintf('Expected Run Time (singleCPU) = %4.1fs, %4.1fm, %4.1fh',t_total,t_total/60, t_total/(60*60)))
p=gcp;Nw=p.NumWorkers;
t_total=t_total/Nw;
disp(sprintf('Expected Run Time (MulCPU) = %4.1fs, %4.1fm, %4.1fh',t_total,t_total/60, t_total/(60*60)))

[t_end_txt,t_left_seconds]=time_loop_end(t0,1,Nd/Nw);
disp(sprintf('Expected end: %s',t_end_txt))

figure(7);clf;
i1=4;i2=6;
subplot(1,2,1);
plot(m_prior{i1},m_prior{i2},'k.');
hold on,
plot(m_post{i1},m_post{i2},'r*');
hold off
xlabel('V_p')
ylabel('\rho')

subplot(1,2,2);
plot(d_prior{1}(1,:),d_prior{1}(2,:),'k.')
hold on,
plot(data{1}.d_obs(1),data{1}.d_obs(2),'g.','MarkerSize',52)
plot(d_post{1}(1,:),d_post{1}(2,:),'r*');
hold off
xlabel('r_o')
ylabel('g')


%% multest prior(post
%id_arr=[1,1000,10000,40000];
rng(1);
id_arr=randomsample(Nd,4);
figure(42);set_paper;
for i=1:length(id_arr);
    id=id_arr(i);
    data=data_mul{id};

    t0=now;
    [logL,evidence,T_est,ABC,dt,iCT]=sippi_abc_logl(ABC,data);
    T=T_est;
    [m_post, P_acc, i_use_all,d_post] = sippi_abc_post_sample(ABC, ns, T, logL);

    data_p=data;
    data_p{1}.Cd=diag([100000,100000,1]);
    %data_p{1}.Cd=diag([1,1,1]/100000);
    [logL_p,evidence_p,T_p]=sippi_abc_logl(ABC,data_p);
    [m_prior, P_acc_prior, i_use_all_prior,d_prior] = sippi_abc_post_sample(ABC, ns,T_p, logL_p);

    subplot(2,ceil(length(id_arr)/2),i);
    plot(d_prior{1}(1,:),d_prior{1}(2,:),'k.','MarkerSize',15)
    hold on,
    plot(data{1}.d_obs(1),data{1}.d_obs(2),'g.','MarkerSize',52)
    plot(d_post{1}(1,:),d_post{1}(2,:),'r.','MarkerSize',12);
    hold off
    xlabel('r_o')
    ylabel('g')

    grid on

    title(sprintf('id=%d, T=%3.1f',id,T))
    drawnow;
end
print_mul('Mina_compare_prior_post_data')

pause(3)

%% MUL INVERSION
Ndata = length(data{1}.d_obs);
D_prior=zeros(Nd,Ndata,Nr);
D_post=zeros(Nd,Ndata,Nr);

R_sat_g=zeros(Nd,Nr);
R_sat_o=zeros(Nd,Nr);
R_sat_b=zeros(Nd,Nr);
R_v_clay=zeros(Nd,Nr);
R_depth=zeros(Nd,Nr);

Rpr_sat_g=zeros(Nd,Nr);
Rpr_sat_o=zeros(Nd,Nr);
Rpr_sat_b=zeros(Nd,Nr);
Rpr_v_clay=zeros(Nd,Nr);
Rpr_depth=zeros(Nd,Nr);

H_sat = zeros(Nd,1);
P_v_clay=zeros(Nd,1);

M_T=zeros(Nd,1);

n_p=101;p = linspace(0,1,n_p);
doComputePriorStat=1;
%parfor id=1:Nd
t0=now;
%doPlot=4;
progress_out(sprintf('%s: Start inversion',mfilename))

parfor id=1:Nd;
    %for id=1:Nd;
    if mod(id,2000)==0,
        [t_end_txt,t_left_seconds]=time_loop_end(t0,id,Nd);
        progress_txt(id,Nd,t_end_txt);
    end
    data = data_mul{id};
    [logL,evidence,T_est]=sippi_abc_logl(ABC,data);
    T=T_est;
    [m_real,P_acc_prior, i_use_all_prior,d_post] = sippi_abc_post_sample(ABC, Nr, T, logL);
    if doComputePriorStat==1;
        data_p=data;
        data_p{1}.Cd=diag([100000,100000,1]);

        [logL_p,evidence_p,T_p]=sippi_abc_logl(ABC,data_p);
        [m_prior, P_acc_prior, i_use_all_prior,d_prior] = sippi_abc_post_sample(ABC, Nr,T_p, logL_p);

        D_prior(id,:,:)=d_prior{1};
        Rpr_sat_g(id,:)=m_prior{1};
        Rpr_sat_o(id,:)=m_prior{2};
        Rpr_sat_b(id,:)=1-(m_prior{2}+m_prior{1});
        Rpr_v_clay(id,:)=m_prior{3};
        Rpr_depth(id,:)=m_prior{7};

        %%
        if (doPlot>3)&(mod(id,10)==0)
            %%
            figure(43);
            subplot(1,2,1)
            plot(D.r0_sim_no_noise(100:100:end), D.g_sim_no_noise(100:100:end),'.','Color',[.8 .8 .8])
            hold on
            plot(d_prior{1}(1,:),d_prior{1}(2,:),'k.','MarkerSize',18)
            scatter(d_prior{1}(1,:),d_prior{1}(2,:),10,m_prior{4},'filled')
            plot(d_post{1}(1,:),d_post{1}(2,:),'r.','MarkerSize',18)
            %scatter(d_post{1}(1,:),d_post{1}(2,:),10,m_real{4},'filled')
            cb=colorbar;set(get(cb,'Ylabel'),'String','C_p')
            plot(data{1}.d_obs(1),data{1}.d_obs(2),'k.','MarkerSize',57)
            plot(data{1}.d_obs(1),data{1}.d_obs(2),'r.','MarkerSize',37)
            plot(data{1}.d_obs(1),data{1}.d_obs(2),'k.','MarkerSize',7)
            hold off
            xlabel('r0');ylabel('g')
            subplot(1,2,2)
            scatter3(m_prior{4},m_prior{5},m_prior{6},10,d_prior{1}(1,:),'filled')
            cb=colorbar;set(get(cb,'Ylabel'),'String','r_0')
            %plot3(m_prior{4},m_prior{5},m_prior{6},'b.')
            hold on
            plot3(m_real{4},m_real{5},m_real{6},'ko')
            hold off
            xlabel('v_p');ylabel('v_s');zlabel('Density')
            box on
            legend('\rho','\sigma')
            title(sprintf('id=%d, T=%3.1f',id, T))
            drawnow;
        end


    end

    M_T(id)=T;

    D_post(id,:,:)=d_post{1};

    R_sat_g(id,:)=m_real{1};
    R_sat_o(id,:)=m_real{2};
    R_sat_b(id,:)=1-(m_real{2}+m_real{1});
    R_v_clay(id,:)=m_real{3};
    R_depth(id,:)=m_real{7};

    %% Useful?
    P_sat = cumsum([R_sat_g(id,:);R_sat_o(id,:);R_sat_b(id,:)]);
    R_cat = zeros(n_p,Nr);
    for ir=1:Nr
        R_cat(:,ir)=1;
        for icat = 1:(size(P_sat,1)-1)
            iic=find(p>P_sat(icat,ir));
            R_cat(iic,ir)=icat+1;
        end
    end

end
progress_out(sprintf('%s: end inversion',mfilename))

M_sat_g=reshape(mean(R_sat_g'),ny,nx);
M_sat_o=reshape(mean(R_sat_o'),ny,nx);
M_sat_b=reshape(mean(1-(R_sat_g'+R_sat_o')),ny,nx);;
M_v_clay=reshape(mean(R_v_clay'),ny,nx);
M_depth=mean(R_depth');
P_v_clay = reshape(sum(R_v_clay<0.1,2)/Nr,ny,nx);
P_oil_50 = reshape(sum(R_sat_o'>0.5)/Nr,ny,nx);
P_gas_50 = reshape(sum(R_sat_g'>0.5)/Nr,ny,nx);

Mpr_sat_g=reshape(mean(Rpr_sat_g'),ny,nx);
Mpr_sat_o=reshape(mean(Rpr_sat_o'),ny,nx);
Mpr_sat_b=reshape(mean(1-(Rpr_sat_g'+Rpr_sat_o')),ny,nx);;
Mpr_v_clay=reshape(mean(Rpr_v_clay'),ny,nx);
Mpr_depth=mean(Rpr_depth');
Ppr_v_clay = reshape(sum(Rpr_v_clay<0.1,2)/Nr,ny,nx);
Ppr_oil_50 = reshape(sum(Rpr_sat_o'>0.5)/Nr,ny,nx);
Ppr_gas_50 = reshape(sum(Rpr_sat_g'>0.5)/Nr,ny,nx);



disp(sprintf('Rejection sampling done'))
d_std = sqrt(diag(data{1}.Cd));
txt = sprintf('Mina_N%d_Nd%d_Nlu%d_Nr%d_%g_%g_%g',N,Nd,Nlu,Nr,1000*d_std(1),1000*d_std(2),d_std(3))
%save(txt,'M_*','P_*','ABC','D*')
save(space2char(txt,'_','\.'),'M_*','P_*', 'R_*','Mpr_*','Ppr_*', 'Rpr_*','D_*')
save([space2char(txt,'_','\.'),'_all'])
%% show results
figure(21);clf
plot(D.r_0_sim, D.g_sim,'.','MarkerSize',.1)
hold on
plot(D.r0_sim_no_noise, D.g_sim_no_noise,'k.','MarkerSize',.1)
plot(D.r_0_data, D.g_data,'r.','MarkerSize',1)
hold off
%%

figure(11);set_paper('landscape')
subplot(2,3,1);
imagesc(x,y,reshape(M_sat_g,ny,nx));axis image;
title('sat_g mean')
caxis([0 1]);colormap(jet)
colorbar

subplot(2,3,2);
imagesc(x,y,reshape(M_sat_o,ny,nx));axis image;
title('sat_o mean')
caxis([0 1]);colormap(jet)
colorbar

subplot(2,3,3);
imagesc(x,y,reshape(1-(M_sat_o+M_sat_g),ny,nx));axis image;
title('sat_b mean')
caxis([0 1]);colormap(jet)
colorbar

subplot(2,3,4);
imagesc(x,y,reshape(M_v_clay,ny,nx));axis image;
title('v_clay mean')
caxis([0 0.5])
colorbar

subplot(2,3,5);
imagesc(x,y,reshape(P_v_clay,ny,nx));axis image;
caxis([0 1]);colormap(gca,flipud(hot))
title('P(v_clay>0.1)')
colorbar


subplot(2,3,6);
imagesc(x,y,reshape(M_T,ny,nx));axis image
colorbar
caxis([1 10]);colormap(gca,flipud(gray))
title('Annealing Temperature')

allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'ydir','normal')

print_mul(txt)

%%
if doComputePriorStat==1
    cmap=cmap_geosoft;
    cmap=[cmap;(copper(2*size(cmap,1)))];
    cmap=cmap_linear([1 1 1; 1 0 0; 1 1 1; 0 0 1; 0 0 0],[0 0.15 0.33, 0.5 1]);
    cmap=cmap_linear([1 1 1; 1 1 1;1 0 0;  0 0 1; 0 0 0],[0 0.05 0.33 0.5 1]);

    figure(12);%set_paper('landscape')

    subplot(3,2,1);
    imagesc(x,y,reshape(Mpr_sat_o,ny,nx));axis image;
    title('sat_o mean - prior')
    caxis([0 1]);
    colormap(gca,cmap)
    colorbar
    subplot(3,2,2);
    imagesc(x,y,reshape(M_sat_o,ny,nx));axis image;
    title('sat_o mean - post')
    caxis([0 1]);
    colormap(gca,cmap)
    colorbar

    subplot(3,2,3);
    imagesc(x,y,reshape(Mpr_sat_g,ny,nx));axis image;
    title('sat_g mean - prior')
    caxis([0 1]);
    colormap(gca,cmap)
    colorbar
    subplot(3,2,4);
    imagesc(x,y,reshape(M_sat_g,ny,nx));axis image;
    title('sat_g mean - post')
    caxis([0 1]);
    colormap(gca,cmap)
    colorbar

    subplot(3,2,5);
    imagesc(x,y,reshape(Mpr_sat_b,ny,nx));axis image;
    title('sat_b mean - prior')
    caxis([0 1]);
    colormap(gca,cmap)
    colorbar
    subplot(3,2,6);
    imagesc(x,y,reshape(M_sat_b,ny,nx));axis image;
    title('sat_b mean - post')
    caxis([0 1]);
    colormap(gca,cmap)
    colorbar

    allAxesInFigure = findall(gcf,'type','axes');
    set(allAxesInFigure,'ydir','normal')

    print_mul([txt,'_satg_prior_post'])

end
%%
figure(41);clf;set_paper;
subplot(2,2,1)
imagesc(x,y,Ppr_oil_50);
title('P(sat_o>0.5), \rho')
axis image;caxis([0,1]);colormap(gca,flipud(hot))
colorbar
subplot(2,2,2)
imagesc(x,y,Ppr_gas_50);
title('P(sat_g>0.5), \rho')
axis image;caxis([0,1]);colormap(gca,flipud(hot))
colorbar

subplot(2,2,3)
imagesc(x,y,P_oil_50);
title('P(sat_o>0.5), \sigma')
axis image;caxis([0,1]);colormap(gca,flipud(hot))
colorbar
subplot(2,2,4)
imagesc(x,y,P_gas_50);
title('P(sat_g>0.5), \sigma')
axis image;caxis([0,1]);colormap(gca,flipud(hot))
colorbar



allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'ydir','normal')

print_mul([txt,'_P_gas_oil'])

%% Information content
for iy=1:ny
    for ix=1:nx
        p_prior=[Mpr_sat_o(iy,ix),Mpr_sat_g(iy,ix),Mpr_sat_b(iy,ix)];
        p_post=[M_sat_o(iy,ix),M_sat_g(iy,ix),M_sat_b(iy,ix)];
        KL(iy,ix)=kl(p_prior,p_post);
        H_prior(iy,ix)=entropy(p_prior,3);
        H_post(iy,ix)=entropy(p_post,3);
    end
end
figure(21);
imagesc(x,y,KL);
colorbar
axis image
colormap(gca,1-gray)
title('KL(prior,post)')
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'ydir','normal')
print_mul([txt,'_kl'])

figure(22);
subplot(1,2,1);
imagesc(x,y,H_prior);
colorbar;axis image;caxis([0 1])
title('H(prior)')
colormap(gca,1-gray)
subplot(1,2,2);
imagesc(x,y,H_post);
colorbar;axis image;caxis([0 1])
title('H(post)')
colormap(gca,1-gray)
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'ydir','normal')
print_mul([txt,'_entropy'])


%%
ii = ceil(linspace(1,N,20000));
figure(31);
plot3(D.g_sim(ii), D.r_0_sim(ii),D.depth_all(ii),'k.','MarkerSize',.1);
hold on
scatter3(D.g_data(:), D.r_0_data(:),D.depth(:),10,M_T(:),'filled');
hold off
caxis([0 10]);
xlabel('g')
ylabel('r_0')
zlabel('depth')
title('T')
colorbar
grid on

print_mul([txt,'_scatter'])

%%
i_p=randomsample(Nd,10000);
figure(32);
subplot(1,2,1)
%plot(M_sat_g(:),M_sat_o(:),'k.')
scatter(M_sat_g(i_p),M_sat_o(i_p),12,M_depth(i_p),'filled')
xlabel('s_g');ylabel('s_o')
axis image;axis([0 1 0 1]);grid on;

subplot(1,2,2)
[Z,x_arr,y_arr] = hist2(M_sat_g(:),M_sat_o(:),21,35);
h=pcolor(x_arr,y_arr,log10(Z'));
shading flat
hold on
i_p = randomsample(Nd,10000);
plot(M_sat_g(i_p),M_sat_o(i_p),'k.','MarkerSize',2)
%scatter(M_sat_g(i_p),M_sat_o(i_p),.1,M_T(i_p),'filled')
hold off
xlabel('s_g');ylabel('s_o')
set(gca,'ydir','normal')
colormap(gca,flipud(hot))
axis image;axis([0 1 0 1]);grid on;
sgtitle('posterior')

%% PRIOR_POST
figure(35);set_paper;

subplot(1,3,1)
p1=plot(D.depth_all(i_p), D.sat_o_log(i_p),'k.');
hold on
p2=plot(M_depth(:),M_sat_o(:),'r.');
hold off
xlabel('Depth');ylabel('S_o')
legend([p1,p2],{'Prior','Post'})

subplot(1,3,2)
p1=plot(D.depth_all(i_p), D.sat_g_log(i_p),'k.');
hold on
p2=plot(M_depth(:),M_sat_g(:),'r.');
hold off
xlabel('Depth');ylabel('S_g')
legend([p1,p2],{'Prior','Post'})

subplot(1,3,3)
p1=plot(D.sat_o_log(i_p), D.sat_g_log(i_p),'k.');
hold on
p2=plot(M_sat_o(:),M_sat_g(:),'r.');
hold off
xlabel('S_o');ylabel('S_g')
legend([p1,p2],{'Prior','Post'})
print_mul([txt,'_priorpost'])
