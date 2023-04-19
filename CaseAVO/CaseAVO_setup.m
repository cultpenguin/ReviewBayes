%%
close all;

if ~exist('useCase','var')
    useCase=3;
end
if ~exist('Cd','var')
    %Cd = [0.0025, -0.005; -.005, 0.0225]; % (r0,g)
    Cd = [0.0025, -0.005 0 ; -.005, 0.0225 0 ; 0 0 1^2]; % (r0,g)
    % Lower the noise, as the noise is probably too high!!!
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.125*sqrt(abs(Cd(1:2,1:2)))).^2;
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.25*sqrt(abs(Cd(1:2,1:2)))).^2;
    %Cd(1:2,1:2) = sign(Cd(1:2,1:2)).*+(0.5*sqrt(abs(Cd(1:2,1:2)))).^2;
end

%%
disp('Load data')
%D=load('data_cube_ALL.mat');
if useCase==1;
    D=load('data_cube_UPDATED.mat'); txt='caseAVO_0';
elseif useCase==2;
    D=load('data_cube_UPDATED2.mat');  txt='caseAVO_2';
else
    D=load('data_cube_UPDATED3_no_corr'); txt='caseAVO_3_no_corr';
end

D.depth_all=repmat(D.depth,[100,1]);
ii=find(D.depth_all>2140 & D.depth_all<2141 );
x=unique(D.xline_obs);nx=length(x);
y=unique(D.inline_obs);ny=length(y);
nxy=nx*ny;

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
subplot(1,3,[1])
scatter(D.xline_obs,D.inline_obs,5,D.depth,'filled')
axis image;colorbar
xlabel('Xline');ylabel('Inline');title('depth')
subplot(1,3,2)
scatter(D.xline_obs,D.inline_obs,5,D.g_data,'filled')
axis image;colorbar
xlabel('Xline');ylabel('Inline');title('g')
subplot(1,3,3)
scatter(D.xline_obs,D.inline_obs,5,D.r_0_data,'filled')
axis image;colorbar
xlabel('Xline');ylabel('Inline');title('R0')
print_mul(sprintf('%s_data',txt))

figure(2);
subplot(2,2,1)
scatter(D.v_p(ii),D.rho(ii),4,D.g_sim(ii),'filled')
xlabel('V_p');ylabel('rho');title('g')
caxis([-1 1].*.3);colormap(gca,cmap_linear)
subplot(2,2,2)
scatter(D.v_p(ii),D.rho(ii),4,D.r_0_sim(ii),'filled')
xlabel('V_p');ylabel('rho');title('R_0')
caxis([-1 1].*.2);colormap(gca,cmap_linear)
colorbar
subplot(2,2,3)
scatter(D.v_p(ii),D.rho(ii),4,D.g_sim_no_noise(ii),'filled')
xlabel('V_p');ylabel('rho');title('g')
caxis([-1 1].*.3);colormap(gca,cmap_linear)
colorbar
subplot(2,2,4)
scatter(D.v_p(ii),D.rho(ii),4,D.r0_sim_no_noise(ii),'filled')
xlabel('V_p');ylabel('rho');title('R_0')
caxis([-1 1].*.2);colormap(gca,cmap_linear)
colorbar

%% plot prior realizations
figure(3);
nr=5;
ns=4;
v_clay_p = zeros(ny,nx);
sat_g_p = zeros(ny,nx);
sat_o_p = zeros(ny,nx);
sat_b_p = zeros(ny,nx);
for i=1:nr;
    i1=(i-1)*nxy+1;
    i2=i*nxy;
    % vclay
    v_clay_p=v_clay_p+reshape(D.v_clay_log(i1:i2),ny,nx);
    subplot(nr,ns,1+(i-1)*ns)
    imagesc(x,y,reshape(D.v_clay_log(i1:i2),ny,nx));
    title(sprintf('v_c, sim #%d',i))
    caxis([0 0.6])
    axis image
    colorbar

    % sat_g
    subplot(nr,ns,2+(i-1)*ns)
    sat_g_p=sat_g_p+reshape(D.sat_g_log(i1:i2),ny,nx);
    imagesc(x,y,reshape(D.sat_g_log(i1:i2),ny,nx));
    caxis([0 1])
    title(sprintf('sat_g, sim #%d',i))   
    axis image

    % sat_o
    subplot(nr,ns,3+(i-1)*ns)
    sat_o_p=sat_g_p+reshape(D.sat_o_log(i1:i2),ny,nx);
    imagesc(x,y,reshape(D.sat_o_log(i1:i2),ny,nx));
    caxis([0 1])
    title(sprintf('sat_o, sim #%d',i))   
    axis image

    % sat_b
    sat_b = reshape(1-D.sat_o_log(i1:i2)-D.sat_g_log(i1:i2),ny,nx);
    subplot(nr,ns,4+(i-1)*ns)
    sat_b_p=sat_b_p+sat_b;
    imagesc(x,y,sat_b);
    caxis([0 1])
    title(sprintf('sat_b, sim #%d',i))   
    axis image
    
end

print_mul(sprintf('%s_prior_reals_org',txt))


%%
figure(4);clf;
subplot(2,2,1);imagesc(x,y,v_clay_p./nr);title('Prior mean, v_c');caxis([0 0.2]);axis image;colorbar
subplot(2,2,2);imagesc(x,y,sat_g_p./nr);title('Sat_g');axis image;caxis([0 1]);colorbar
subplot(2,2,3);imagesc(x,y,sat_o_p./nr);title('Sat_o');axis image;caxis([0 1]);colorbar
subplot(2,2,4);imagesc(x,y,sat_b_p./nr);title('Sat_b');axis image;caxis([0 1]);colorbar
print_mul(sprintf('%s_prior_mean_org',txt))

