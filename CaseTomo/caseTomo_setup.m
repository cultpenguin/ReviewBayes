close all
if ~exist('dx','var')
    %dx=0.5;
    dx=0.25;
    dx=0.1;
end
forward.null=[];

ax=[-1 6 0 13];

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

% bimodal distribution
Nd=10000;
prob_chan=0.5;
dd=.015;
d1=randn(1,ceil(Nd*(1-prob_chan)))*.0045+0.145-dd;  %0.1125;
d2=randn(1,ceil(Nd*(prob_chan)))*.0045+0.145+dd; %0.155;
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
prior{im}.cmap=cmap;
prior{im}.o_nscore=o_nscore;

%% SETUP THE FORWARD MODEL
forward.linear_m=m0; % Needed when m0 is set to 0 (using target distribution, and forward.linear=1;)�
forward.sources=D.S;
forward.receivers=D.R;
if ~isfield(forward,'type');
    forward.type='eikonal';
    forward.type='fat';
end
if ~isfield(forward,'linear'); forward.linear=1;end
if ~isfield(forward,'freq'); forward.freq=0.1;end
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
    Nme=600;
    [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward,prior,Nme);
    
    % ASSIGN MODELING ERROR TO DATA
    for id=1:length(data);
        data{id}.dt=dt{id};
        data{id}.Ct=Ct{id};
    end
else
    forward_full.type='none';
end

% 

txt=sprintf('caseBayesian_dx%d_F%s-%s_ME%d',ceil(100*dx),forward.type,forward_full.type,comp_model_error);
disp(txt)
save(txt)

%%
figure(1)
for i=1:5;
    subplot(1,5,i);
    m=sippi_prior(prior);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis(cax);colormap(cmap)
    axis image
    axis(ax)
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

%% 
figure(10);
subplot(1,5,1)
for i=1:size(D.S,1)
    plot([D.S(i,1),D.R(i,1)],[D.S(i,2),D.R(i,2)],'ko-','LineWidth',.1,'MarkerSize',1)
    hold on
end
hold off
axis image
axis(ax)
xlabel('X (m)');ylabel('Y (m)')
print -dpng -r300 caseTomo_SR

%% Frechet
dv=0.0001;
iray=10;
m_ref=sippi_prior(prior);
d_ref=sippi_forward(m_ref,forward,prior);
for ix = 1:length(prior{1}.x);
    progress_txt(ix,length(prior{1}.x),'X')
for iy = 1:length(prior{1}.y);
    m_test=m_ref;
    m_test{1}(iy,ix)=m_ref{1}(iy,ix)+dv;
    d_test=sippi_forward(m_test,forward,prior);
    m_fre(iy,ix)= (d_ref{1}(iray)-d_test{1}(iray))/dv;
end
end

%%
figure(21);clf
subplot(1,2,1);
imagesc(prior{1}.x, prior{1}.y, m_ref{1});
axis image;
caxis(prior{1}.cax);colormap(gca,prior{1}.cmap)
title('Example velocity model')
subplot(1,2,2);
imagesc(prior{1}.x, prior{1}.y, m_fre);
colormap(gca,cmap_linear)
cax=caxis;
caxis([-1 1].*max(abs(cax)))
xlabel('X (m)');ylabel('Y (m)');
hold on;
plot(D.S(iray,1),D.S(iray,2),'r*','MarkerSize',10);
plot(D.R(iray,1),D.R(iray,2),'ro','MarkerSize',10)
hold off
axis image;
colorbar_shift;
title(sprintf('Frechét derivative [dd/dm]'))
print('-dpng','-r300',sprintf('caseTomo_frechet_%s_L%d_SR%d',forward.type,forward.linear,iray))




    


