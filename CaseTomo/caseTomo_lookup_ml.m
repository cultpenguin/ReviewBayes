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
if ~isfield(prior{1},'cax');prior{1}.cax=[0.1150    0.1750];end
if ~isfield(prior{1},'cax_std');prior{1}.cax_std=[0 0.02];end
if ~isfield(prior{1},'cmap');prior{1}.cmap=jet;end


txt_out = sprintf('%s_lu_noise_N%d',txt,N);
disp(txt_out)


if exist([txt_out,'.mat'],'file')
    disp(sprintf('Loading %s',txt_out))
    load(txt_out,'ABC')
else
    ABC.useParfor=1;
    ABC.simulate_noise=1;ABC.data=data; % We need to simulate NOISE
    
    ABC=sippi_abc_setup(prior,forward,N,Nme,ABC);
    save(txt_out,'ABC','-v7.3')
end


% get m{2}:AREA and m{3};P(v<vmax)
dx=prior{1}.x(2)-prior{1}.x(1);
for i=1:N
    if forward.is_slowness==1;
        v=1./ABC.m{i}{1}(:)';
    else
        v=ABC.m{i}{1}(:)';
    end
    [P,A,vmax]=caseTomo_reals_to_P_area(v,dx);
    
    ABC.m{i}{2}=A;
    ABC.m{i}{3}=ABC.m{i}{1}<vmax;
end

%% Export to hdf5
[h5file,M,D,Dobs]=sippi_abc_to_h5(ABC,txt_out);


%% VELOCITY
clear ml
ml.type = 'regression';
ml.hidden_layers = 6;
ml.hidden_units = 40;
ml.use_log = 0;
ml.MaxEpochs=500;
ml.MiniBatchSize=4*128;
ml.MaxIteNotImproving=10;
%ml.ExecutionEnvironment = 'cpu';
ml.normalize=1; % normalize data
% ml.id=1; % noise free data
ml.id=2; % noisy data
ml.im=1; % velocity
ml.Plots = 'training-progress';
ml.Plots = 'none';

[net,ml,ABC]=sippi_abc_ml_setup(ABC,ml);
%%
D=data{1}.d_obs;
%D=ABC.d{3}{2}
M_est=sippi_abc_ml_predict(ABC,D);
x=prior{1}.x;
y=prior{1}.y;
m_est=reshape(M_est,length(y),length(x));

figure(41);clf
subplot(1,3,1)
imagesc(x,y,m_est)
axis image;colormap(prior{1}.cmap)
caxis(prior{1}.cax)
title('\sigma(m) - mean')
colorbar
print_mul(sprintf('%s_post_mean',txt_out))


%% AREA
iml=2;
ml.type = 'regression';
ml.hidden_layers = 4;
ml.hidden_units = 100;
ml.use_log = 0;
ml.MiniBatchSize=4*128;
ml.MaxIteNotImproving=10;
ml.normalize=1; % normalize data
ml.id=2; % noise data
ml.im=iml; % AREA
[net,ml,ABC]=sippi_abc_ml_setup(ABC,ml);
D=data{1}.d_obs;
A_est=sippi_abc_ml_predict(ABC,D,iml)

%% TEST
for i=1:1000;
    D=ABC.d{i}{2};
    A_est=sippi_abc_ml_predict(ABC,D,iml);
    [A_est,ABC.m{i}{2}];
    A_est_all(i) = A_est;
    A_true_all(i) = ABC.m{i}{2};
end
plot(A_est_all,A_true_all,'.')
xlabel('A_{est}')
ylabel('A_{true}')



%% ABC
[logL,evidence,T_est]=sippi_abc_logl(ABC,data);
%%
ns=100;
T=8;T_est;
[i_use_all,P_acc]=sippi_abc_post_sample_logl(logL,ns,T);
clear m1 m2 m3
for i=1:ns;
    m1_prior(:,:,i)=ABC.m{i_use_all(i)}{1};
    m2_prior(i)=ABC.m{i_use_all(i)}{2};
    m3_prior(:,:,i)=ABC.m{i_use_all(i)}{3};
    m1_post(:,:,i)=ABC.m{i_use_all(i)}{1};
    m2_post(i)=ABC.m{i_use_all(i)}{2};
    m3_post(:,:,i)=ABC.m{i_use_all(i)}{3};
end
[em,ev,emode]=etype(m1_post);
[P]=etype(m3_post);
figure(21);clf;
subplot(1,3,1);
imagesc(x,y,em);axis image
caxis(prior{1}.cax)
colorbar
subplot(1,3,2);
imagesc(x,y,sqrt(ev));axis image
colorbar
subplot(1,3,3);
imagesc(x,y,P);axis image
colorbar
nx=length(prior{1}.x);
ny=length(prior{1}.y);
nxy=nx*ny;
txt_out2 = sprintf('%s_T%d',txt_out,ceil(T));
caseTomo_plot_post_stats(reshape(m1_post,[nxy,ns])',reshape(m1_prior,[nxy,ns])',prior,txt_out2);
%%
return
%% PROBABILITY
iml=3;
ml.type = 'classification';
ml.hidden_layers = 2;
ml.hidden_units = 100;
ml.use_log = 0;
ml.MiniBatchSize=4*128;
ml.MaxIteNotImproving=10;
ml.normalize=1; % normalize data
ml.id=2; % noise data
ml.im=iml; % 
[net,ml,ABC]=sippi_abc_ml_setup(ABC,ml);
D=data{1}.d_obs;
M_est=sippi_abc_ml_predict(ABC,D,iml);
