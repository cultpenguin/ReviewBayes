close all
if ~exist('rseed','var');rseed=1;end
if ~exist('dx','var')
    %dx=0.5;
    dx=0.25;
    %dx=0.1;
end
if ~exist('doComputeFrechet','var')
    doComputeFrechet=0;
end
if ~exist('useCase','var')
    useCase='Kallerup';
    %useCase='Arrenaes';
end
%%
if rseed>0; rng('default');rng(rseed);end
forward.null=[];

cmap=jet;

% Load data
if strcmp(useCase,'Kallerup')
    addpath Kallerup
    %K=load('KallerupJensenOutput');
    K=load('KallerupJensenOutput_small','ant_pos','data','prior','c0');
    options.txt='Kallerup';
    ax=[-.5 4.0 0 7];
    cax=[0.07 0.18];
    D.S=K.ant_pos(:,1:2);
    D.R=K.ant_pos(:,3:4);
    D.d_obs = K.data{1}.d_obs;
    D.d_std = K.data{1}.d_std;
    D.Ct = K.data{1}.Ct;
    D.dt = K.data{1}.dt;
    D.Ct = K.data{1}.CD;
    addExtra=1;
    if addExtra==1
        D.Ct=D.Ct+4;
    end

    %D.dt = K.data{1}.dD;
    clear data
else
    D=load('AM13_data.mat');
    options.txt='AM13_bimodal';
    ax=[-1 6 0 13];
    cax=.145+[-1 1].*.03;

end




%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=D.Ct; % modelization and static error
if isfield(D,'dt');
    data{id}.dt=D.dt; % modelization and static error
end



%% SETUP PRIOR

if strcmp(useCase,'Kallerup')

    p = K.prior;
    clear prior

    % Use velocity as prior
    d_target_eps = p{1}.d_target;
    d_target_vel= (sqrt(K.c0.^2./d_target_eps)*10^-9);
    d=d_target_vel;
    n=9000;

    useAltHigh=0;
    if useAltHigh==1;
        v1=d(d<0.1);
        v2=d(d>0.1);
        v2a=randn(n,1)*1.5*std(v2)+0.18;
        d=[v1;v2a];
    end

    var0=var(d);
    m0=mean(d);
    %d=randn(2*n,1)*sqrt(var0)+m0;d(d<0.02)=0.02;

    [d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
    im=1;
    clear prior
    prior{im}.d_target=d;
    prior{im}.o_nscore=o_nscore;
    prior{im}.type='FFTMA';
    prior{im}.name='Velocity (m/ns)';
    prior{im}.Va=sprintf(sprintf('%5.4f Sph(15,90,.1)',var0));
    %prior{im}.Va=sprintf(sprintf('%5.4f Sph(5,90,.3)',var0));
    prior{im}.x=[ax(1):dx:ax(2)];
    prior{im}.y=[ax(3):dx:ax(4)];
    prior{im}.cax=cax;
    prior{im}.cmap=cmap;

else

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
    prior{im}.d_target=d;
    prior{im}.type='FFTMA';
    prior{im}.name='Velocity (m/ns)';
    prior{im}.Va=sprintf(sprintf('%5.4f Sph(6,90,.2)',var0));
    %prior{im}.x=[-1:dx:6];
    %prior{im}.y=[0:dx:13];
    prior{im}.x=[ax(1):dx:ax(2)];
    prior{im}.y=[ax(3):dx:ax(4)];

    prior{im}.cax=cax;
    prior{im}.cmap=cmap;
    prior{im}.o_nscore=o_nscore;

end

%% SETUP THE FORWARD MODEL
forward.linear_m=m0; % Needed when m0 is set to 0 (using target distribution, and forward.linear=1;)�
forward.sources=D.S;
forward.receivers=D.R;
if ~isfield(forward,'type');
    forward.type='eikonal';
    forward.type='fat';
    forward.type='ray';
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

txt=sprintf('caseTomo_%s_dx%d_F%s-%s_ME%d',useCase,ceil(100*dx),forward.type,forward_full.type,comp_model_error);
disp(txt)
save(txt)

%%
figure(1)
if rseed>0; rng('default');rng(rseed);end
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
cb=colorbar_shift;
set(get(cb,'Ylabel'),'String','Velocity (m/ns)')
print_mul(sprintf('%s_%s',txt,'priorsample'))

try
    figure;
    histogram(prior{1}.o_nscore.d,linspace(cax(1),cax(2),31),'Normalization','Probability');
    xlabel('Velocity (m/\mus)')
    ylabel('Probability')
    grid on
    print_mul(sprintf('%s_%s',txt,'prior1Dmarg'))
end

%%
%m=sippi_prior(prior);
d=sippi_forward(m,forward,prior,data);
figure(8);clf;
%subplot(1,2,1)
%imagesc(prior{1}.x,prior{1}.y,m{1});
%caxis(cax);colormap(cmap)
%axis image
%axis(ax)
%xlabel('X (m)')
%ylabel('Y (m)')
%
%subplot(1,2,2)
plot(data{1}.d_obs,'k-','LineWidth',1.5)
hold on
plot(d{1},'r-','LineWidth',1)
hold off
xlabel('Data #')
ylabel('travel time (ns)')
grid on
legend({'d_{obs}','g(m)'})
print_mul(sprintf('%s_%s',txt,'_forward'))



%% 'ray' Coverage
for i=1:size(D.S,1)
    dis(i)=edist(D.S(i,:),D.R(i,:));
    vapp(i)=dis(i)./data{1}.d_obs(i);
end
figure(9);clf
plot(D.S(:,1),D.S(:,2),'k.')
hold on
plot(D.R(:,1),D.R(:,2),'k.')
colormap(gca,cmap);caxis(cax);
cb=colorbar;
set(get(cb,'Ylabel'),'String','Velocity (m/ns)')


for i=1:size(D.S,1)
    %icol=ceil(interp1([cax],[1 size(cmap,1)],vapp(i),'nearest'))
    icol(i)=ceil(interp1([cax],[1 size(cmap,1)],vapp(i),'linear','extrap'));
    if vapp(i)<cax(1);icol(i)=1;end
    if vapp(i)>cax(2);icol(i)=size(cmap,1);end
    lw=2*(vapp-0.1)./(0.1);
    lw(lw<0.01)=0.01;
    plot([D.S(i,1),D.R(i,1)],[D.S(i,2),D.R(i,2)],'-','LineWidth',lw(i),'MarkerSize',1,'Color',cmap(icol(i),:))
    %plot([D.S(i,1),D.R(i,1)],[D.S(i,2),D.R(i,2)],'-','LineWidth',.1)
end
grid on
hold off
axis image
axis(ax)
set(gca,'ydir','reverse')
print_mul(sprintf('%s_%s',txt,'vapp'))





%% Frechet
if doComputeFrechet==1;

    dv=0.0005;
    iray=10;
    if rseed>0; rng(rseed);end
    m_ref=sippi_prior(prior);
    d_ref=sippi_forward(m_ref,forward,prior);
    m_fre=zeros(length(prior{1}.y),length(prior{1}.x))
    for ix = 1:length(prior{1}.x);
        progress_txt(ix,length(prior{1}.x),'X')
        parfor iy = 1:length(prior{1}.y);
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
    print_mul(sprintf('%s_frechet_%s_L%d_SR%d',txt,forward.type,forward.linear,iray))
end






