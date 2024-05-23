%%
%clear all;close all;
doPlot=1;
%wx=0; % with of local Nhood X-dur
%wy=0; % width of local Nhood Y-dir
%dx=1;
%dy=1;


if ~exist('Nlu','var')
    Nlu=1000000;
    Nlu=500000;
    Nlu=100000;
    Nlu=1000;
    %Nlu=100;
end
if ~exist('Nr','var')
    Nr=400;
end
if ~exist('useUncond','var')
    useUncond=1;
end

if ~exist('Cd','var')
    %Cd = [0.0025, -0.005; -.005, 0.0225]; % (r0,g)
    Cd = [0.0025, -0.005 0 ; -.005, 0.0225 0 ; 0 0 1^2]; % (r0,g)
    % Lower the noise, as the noise is probably too high!!!
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.125*sqrt(abs(Cd(1:2,1:2)))).^2;
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.25*sqrt(abs(Cd(1:2,1:2)))).^2;
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.5*sqrt(abs(Cd(1:2,1:2)))).^2;
end

if ~exist('dx','var');dx=1;end
if ~exist('dy','var');dy=dx;end
if ~exist('wx','var');wx=0;end
if ~exist('wy','var');wy=wx;end



%
disp('Load data')
%D=load('data_cube_ALL.mat');
%D=load('data_cube_UPDATED.mat');  Nlu=Nlu+1;
%D=load('data_cube_UPDATED2.mat');    Nlu=Nlu+2; txt = 'no_well_cond';
if useUncond == 1
    progress_out('--> load Uncond Data')
    D=load('data_cube_UNCOND_1000_realisations.mat');D.v_p = D.sat_g;D.v_s =D.sat_g;D.rho = D.sat_g;Nlu=Nlu+3; txt= 'well_uncond'
else
    progress_out('--> load Data')
    D=load('data_cube_1000_realisations.mat');D.v_p = D.sat_g;D.v_s =D.sat_g;D.rho = D.sat_g;Nlu=Nlu+3; txt= 'well_cond'
end

Nd = length(D.g_data);
Nr_prior=length(D.g_sim)/length(D.g_data);
x=unique(D.xline_obs);nx=length(x);
y=unique(D.inline_obs);ny=length(y);
nxy=nx*ny;

txt = sprintf('Mina_wx%d_wy%d_Nr%d_%s',wx,wy,Nr_prior,txt);

%%
%    R_sat_g(id,:)=m_real{1};
%    R_sat_o(id,:)=m_real{2};
%    R_sat_b(id,:)=1-(m_real{2}+m_real{1});
%    R_v_clay(id,:)=m_real{3};

% M
v_clay_log = reshape(D.v_clay_log,[ny nx Nr_prior]);
sat_g_log = reshape(D.sat_g_log,[ny nx Nr_prior]);
sat_o_log = reshape(D.sat_o_log,[ny nx Nr_prior]);

% D
r_0_data = reshape(D.r_0_data,[ny nx]);
g_data = reshape(D.g_data,[ny nx]);
r0_sim_no_noise = reshape(D.r0_sim_no_noise,[ny nx Nr_prior]);
g_sim_no_noise = reshape(D.g_sim_no_noise,[ny nx Nr_prior]);

M_sat_o_prior = mean(sat_o_log,3);%mean(reshape(D.sat_o_log,[Nd,Nr_prior])');
M_sat_g_prior = mean(sat_g_log,3);%mean(reshape(D.sat_g_log,[Nd,Nr_prior])');
M_v_clay_prior = mean(v_clay_log,3);%mean(reshape(D.v_clay_log,[Nd,Nr_prior])');

figure(11);clf;
for i=1:min([1  Nr_prior]);
    i1=(i-1)*nx*ny+1;i2=i*nx*ny;
    subplot(2,3,1);imagesc(v_clay_log(:,:,i));caxis([0 1]*1);axis image;;title('v_c')
    subplot(2,3,2);imagesc(sat_g_log(:,:,i));caxis([0 1]*1);axis image;title('sat_g')
    subplot(2,3,3);imagesc(sat_o_log(:,:,i));caxis([0 1]*1);axis image;title('sat_o')

    subplot(2,3,4);imagesc(reshape(M_v_clay_prior,ny,nx));caxis([0 1]*1);axis image;
    subplot(2,3,5);imagesc(reshape(M_sat_g_prior,ny,nx));caxis([0 1]*1);axis image;
    subplot(2,3,6);imagesc(reshape(M_sat_o_prior,ny,nx));caxis([0 1]*1);axis image;
    drawnow;
end
sgtitle(txt,'interpreter','none')
print_mul([txt,'_prior'])



%% Setup ABC
progress_out('--> SetupABC')
clear ABC

R_sat_g=zeros(Nd,Nr).*NaN;
R_sat_o=zeros(Nd,Nr).*NaN;
R_sat_b=zeros(Nd,Nr).*NaN;
R_v_clay=zeros(Nd,Nr).*NaN;
R_depth=zeros(Nd,Nr).*NaN;

H_sat = zeros(Nd,1).*NaN;
P_v_clay=zeros(Nd,1).*NaN;

M_T=zeros(Nd,1).*NaN;
Nabove=zeros(Nd,1).*NaN;

Nlu = min([Nlu,Nr_prior]);
N=length(D.g_sim);

%%

progress_out(sprintf('--> Start local sampling, wx=%d',wx))

useAnneal=1;
t0=now;
id=0;
for ix=(1+wx):dx:(nx-wx);
    if (ix>(1+wx))&&(mod(ix,1)==0)
        [t_end_txt,t_left_seconds]=time_loop_end(t0,id,Nd);
        progress_txt(id,Nd,t_end_txt);
        figure_focus(21);clf
        subplot(1,2,1)
        imagesc(reshape(mean(R_sat_g'),ny,nx));axis image
        colormap(gca,jet)
        subplot(1,2,2)
        imagesc(reshape(M_T,ny,nx));axis image
        caxis([0.5 5])
        colormap(gca,hot)
        drawnow
    end

    for iy=(1+wy):dy:(ny-wy);
        id=id+1;
        id=(ix-1)*ny+iy;
        clear ABC i k data
        i_use = 1:Nlu;
        iix = [ (ix-wx):1:(ix+wx) ];
        iiy = [ (iy-wy):1:(iy+wy) ];
        [ixx,iyy]=meshgrid(iix,iiy);
        n_local = length(iix)*length(iiy);
        i0 = ceil(n_local/2);
        i=0;
        for k=1:Nlu;
            i=i+1;
            % m
            j=0;

            % HER BEHØVES JO KUN DEN CENTRALE CELLE
            j=j+1;m1=sat_g_log(iiy,iix,k);ABC.m{i}{j}=m1(:);
            j=j+1;m2=sat_o_log(iiy,iix,k);ABC.m{i}{j}=m2(:);
            j=j+1;m3=v_clay_log(iiy,iix,k);ABC.m{i}{j}=m3(:);

            if i==1;
                ABC.m=repmat(ABC.m,[length(i_use),1]);
            end

            % d
            %j=0;
            %j=j+1;ABC.d{i}{j}=[ r0_sim_no_noise(iy,ix,k), g_sim_no_noise(iy,ix,k)]';
            for j=1:length(ixx(:));
                ABC.d{i}{j}=[ r0_sim_no_noise(iyy(j),ixx(j),k), g_sim_no_noise(iyy(j),ixx(j),k)]';
            end
            if i==1;
                ABC.d=repmat(ABC.d,[length(i_use),1]);
            end

        end

        Nlu=length(ABC.m);
        for i=1:n_local
            %data{i}.d_obs = [D.r_0_data(id), D.g_data(id)]';
            data{i}.d_obs = [r_0_data(iyy(i),ixx(i)), g_data(iyy(i),ixx(i))]';
            %data{i}.d_std=sqrt(diag(Cd));
            data{i}.d_std = [0.0000001,0.000000001]';
            data{i}.Cd = Cd(1:2,1:2);
        end
        %%
        % NEXT LINE: MAKE SURE ALL DATA ARE USED!!!!
        ABC.use_sippi_likelihood=1;
        [logL,evidence,T_est]=sippi_abc_logl(ABC,data);
        Nabove(id)=sum((exp(logL-max(logL))>0.1));
        if useAnneal==1;
            T=T_est;
	else
	    T=1;
        end	    
        [m_real,P_acc_prior, i_use_all_prior,d_post] = sippi_abc_post_sample(ABC, Nr, T, logL);

        % POST STAT
        M_T(id)=T;

        %D_post(id,:,:)=d_post{1};

        R_sat_g(id,:)=m_real{1}(i0,:);
        R_sat_o(id,:)=m_real{2}(i0,:);
        R_sat_b(id,:)=1-(m_real{2}(i0,:)+m_real{1}(i0,:));
        R_v_clay(id,:)=m_real{3}(i0,:);



    end
end
progress_out(sprintf('--> End local sampling, wx=%d',wx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE
disp(sprintf('Rejection sampling done'))
d_std = sqrt(diag(data{1}.Cd));
txt = sprintf('Mina_N%d_Nd%d_Nlu%d_Nr%d_wx%d_wy%d_%g_%g_A%d_UC%d',N,Nd,Nlu,Nr,wx,wy,1000*d_std(1),1000*d_std(2),useAnneal,useUncond)
%save(txt,'M_*','P_*','ABC','D*')
  save(space2char(txt,'_','\.'),'M_*','P_*', 'R_*','N*')
%save([space2char(txt,'_','\.'),'_all'])

%% plot
M_sat_g=reshape(mean(R_sat_g'),ny,nx);
M_sat_o=reshape(mean(R_sat_o'),ny,nx);
M_sat_b=reshape(mean(1-(R_sat_g'+R_sat_o')),ny,nx);;
M_v_clay=reshape(mean(R_v_clay'),ny,nx);
M_depth=mean(R_depth');
P_v_clay = reshape(sum(R_v_clay<0.1,2)/Nr,ny,nx);
P_oil_50 = reshape(sum(R_sat_o'>0.5)/Nr,ny,nx);
P_gas_50 = reshape(sum(R_sat_g'>0.5)/Nr,ny,nx);


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




