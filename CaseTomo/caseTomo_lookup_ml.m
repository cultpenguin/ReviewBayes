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
        ABC.useParfor=1;
        ABC.simulate_noise=1;ABC.data=data; % We need to simulate NOISE
        
        ABC=sippi_abc_setup(prior,forward,N,Nme,ABC);
        save(txt_out_nn,'ABC')
    end

end

% get m{2}:AREA and m{3};P(v<vmax)
dx=prior{1}.x(2)-prior{1}.x(1);
for i=1:N
    [P,A,vmax]=caseTomo_reals_to_P_area(ABC.m{1}{1}(:)',dx);
    
    ABC.m{i}{2}=A;
    ABC.m{i}{3}=ABC.m{i}{1}<vmax;
end

%% Export to hdf5
[h5file,M,D,Dobs]=sippi_abc_to_h5(ABC,txt_out);


%% VELOCITY
ml.type = 'regression';
ml.hidden_layers = 2;
ml.hidden_units = 100;
ml.use_log = 0;
ml.MiniBatchSize=4*128;
ml.MaxIteNotImproving=10;
ml.normalize=1; % normalize data
ml.id=2; % noise data
ml.im=1; % velocity
[net,ml,ABC]=sippi_abc_ml_setup(ABC,ml);
%%
D=data{1}.d_obs;
M_est=sippi_abc_ml_predict(ABC,D);
x=prior{1}.x;
y=prior{1}.y;
m_est=reshape(M_est,length(y),length(x));

figure(41);clf
subplot(1,3,1)
imagesc(x,y,m_est)
colormap(prior{1}.cmap)
caxis(prior{1}.cax)
axis image
title('\sigma(m) - mean')
print_mul(sprintf('%s_post_mean',txt_out))


%% AREA
iml=2;
ml.type = 'regression';
ml.hidden_layers = 2;
ml.hidden_units = 100;
ml.use_log = 0;
ml.MiniBatchSize=4*128;
ml.MaxIteNotImproving=10;
ml.normalize=1; % normalize data
ml.id=2; % noise data
ml.im=iml; % AREA
[net,ml,ABC]=sippi_abc_ml_setup(ABC,ml);
D=data{1}.d_obs;
A_est=sippi_abc_ml_predict(ABC,D,iml);


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
