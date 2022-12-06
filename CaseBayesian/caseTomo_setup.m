clear all;close all

cax=0.145+[-1 1].*.03;
cmap=jet;

% Load data
D=load('AM13_data.mat');
options.txt='AM13_bimodal';

%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=D.Ct; % modelization and static error


%% SETUP PRIOR
dx=0.25;
%dx=0.5;
%dx=0.1;

% bimodal distribution
N=10000;
prob_chan=0.5;
dd=.015;
d1=randn(1,ceil(N*(1-prob_chan)))*.0045+0.145-dd;  %0.1125;
d2=randn(1,ceil(N*(prob_chan)))*.0045+0.145+dd; %0.155;
d=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
var0=var(d);
m0=mean(d);

im=1;
clear prior
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.Va=sprintf(sprintf('%5.4f Sph(6,90,.2)',var0));
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=cax;
prior{im}.o_nscore=o_nscore;

m=sippi_prior(prior);
imagesc(m{1});colorbar

%% SETUP THE FORWARD MODEL
forward.linear_m=m0; % Needed when m0 is set to 0 (using target distribution, and forward.linear=1;)�
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
%forward.type='ray_2d';forward.r=2;
forward.type='fat';
forward.linear=1;
forward.freq=0.1;
forward.forward_function='sippi_forward_traveltime';

m=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward,prior,data);

comp_model_error=0;
if comp_model_error==1;
    
    % SETUP THE 'OPTIMAL' FORWARD MODEL
    forward_full.forward_function='sippi_forward_traveltime';
    forward_full.sources=D.S;
    forward_full.receivers=D.R;
    forward_full.type='fat';forward_full.linear=0;forward_full.freq=0.1;
    
    % COMPUTE MODELING ERROR DUE TO USE OF forward AS OPPOSED TO forward_full
    N=600;
    [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward,prior,N);
    
    % ASSIGN MODELING ERROR TO DATA
    for id=1:length(data);
        data{id}.dt=dt{id};
        data{id}.Ct=Ct{id};
    end
else
    forward_full.type='none';
end

% 

txt=sprintf('caseBayesian_dx%d_F%s-%s_ME%d',ceil(100*dx),forward.type,forward_full.type,comp_model_error)
save(txt)

%%
figure(1)
for i=1:5;
    subplot(1,5,i);
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis(cax);colormap(cmap)
    axis image
    xlabel('X (m)')
    if i==1,    ylabel('Y (m)');end
end
print -dpng -r300 caseTomo_priorsample

try
    figure;
    histogram(prior{1}.o_nscore.d)
    xlabel('Velocity (m/\mus)')
    print -dpng -r300 caseTomo_prior1Dmarg
end