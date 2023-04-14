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

%% Setup data
disp('Setup data_mul')
Nd = length(D.g_data);
for i=1:Nd;
    data{1}.d_obs = [D.r_0_data(i), D.g_data(1),D.depth(i)]';
    % NEXT LINE SETS THE EXPECTED STD OF THE DATA NOISE!
    %data{1}.d_std = [0.005,0.005, 1]';
    data{1}.d_std=sqrt(diag(Cd));
    data{1}.d_std = [0.0000001,0.000000001, .00000001]';
    data{1}.Cd = Cd;

    data_mul{i}=data;
end

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