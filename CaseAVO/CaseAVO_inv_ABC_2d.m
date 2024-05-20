%clear
% all;close all;
%useCase=3;Nlu=100000;Nr=400;useRejection=0;
%useCase=2;Nlu=100000;Nr=400;useRejection=1;
%useCase=4;Nlu=100000;Nr=400;useRejection=1;
%useCase=3;Nlu=100000;Nr=400;useRejection=1;

%%
if ~exist('useCase','var');useCase=3;end
if ~exist('Nlu','var');Nlu=1001;end
if ~exist('Nr','var');Nr=400;end
if ~exist('useRejection','var');useRejection=1;end

%% Setup ABC
progress_out('--> START CaseAVO_setup')
CaseAVO_setup
progress_out('--> END CaseAVO_setup')


%%
Nlu=nr;CaseAVO_setup_ABC;
D_r = reshape(D.r0_sim_no_noise,[nx*ny,100]);
D_g = reshape(D.g_sim_no_noise,[nx*ny,100]);

D_sat_g_log = reshape(D.sat_g_log,[nx*ny,100]);
D_sat_o_log = reshape(D.sat_o_log,[nx*ny,100]);
D_v_clay_log = reshape(D.v_clay_log,[nx*ny,100]);

%Nr=10;
R_sat_g=zeros(Nd,Nr);
R_sat_o=zeros(Nd,Nr);
R_sat_b=zeros(Nd,Nr);
R_v_clay=zeros(Nd,Nr);

%%
time_sampling_0=now;
progress_out('--> START SAMPLING')
for id=1:length(data_mul)
    if mod(id,1000)==0
        progress_txt(id,Nd);
    end
    
    data = data_mul{1};
    % [r,g]
    i_use = 1:2;
    data{1}.d_obs = data{1}.d_obs(i_use);
    data{1}.d_std = data{1}.d_std(i_use);
    data{1}.Cd = data{1}.Cd(i_use,i_use);
    
    % setup ABC
    %ABC.d{id{1}}
    for ir=1:nr;
        ABC.d{ir}{1}=[D_r(id,ir), D_g(id,ir)]';

        j=0;
        j=j+1;ABC.m{ir}{j}=D_sat_g_log(id,ir);
        j=j+1;ABC.m{ir}{j}=D_sat_o_log(id,ir);
        j=j+1;ABC.m{ir}{j}=D_v_clay_log(id,ir);
        %j=j+1;ABC.m{ir}{j}=D.v_p(k);
        %j=j+1;ABC.m{ir}{j}=D.v_s(k);
        %j=j+1;ABC.m{ir}{j}=D.rho(k);
        %j=j+1;ABC.m{ir}{j}=D.depth_all(k);
        %j=j+1;ABC.m{i}{j}=D.porosity(k);
    
    end

    %[logL,evidence,T_est,ABC,dt,iCT]=sippi_abc_logl(ABC,data,N);
    [logL]=sippi_abc_logl(ABC,data);
    T=.1;
    ns=Nr;
    [i_use_all,P_acc]=sippi_abc_post_sample_logl(logL,ns,T);


    plData=0;
    if plData==1;
        D1 = [D_r(id,:); D_g(id,:)]';
        figure(1);
        plot(D1(:,1),D1(:,2),'k.','MarkerSize',3);
        hold on;
        plot(D1(i_use_all,1),D1(i_use_all,2),'r.');
        plot(data{1}.d_obs(1),data{1}.d_obs(2),'r*');
        hold off
        drawnow;
    end

    R_sat_g(id,:)=D_sat_g_log(id,i_use_all);
    R_sat_o(id,:)=D_sat_o_log(id,i_use_all);
    R_sat_b(id,:)=1-(D_sat_g_log(id,i_use_all)+D_sat_o_log(id,i_use_all));
    R_v_clay(id,:)=D_v_clay_log(id,i_use_all);

end
progress_out('--> END SAMPLING')
time_sampling_1=now;
time_sampling=(time_sampling_1-time_sampling_0)*3600*24;
%%
M_sat_g=reshape(mean(R_sat_g'),ny,nx);
M_sat_o=reshape(mean(R_sat_o'),ny,nx);
M_sat_b=reshape(mean(1-(R_sat_g'+R_sat_o')),ny,nx);;
M_v_clay=reshape(mean(R_v_clay'),ny,nx);


disp(sprintf('Sampling done in %5.1 minites',time_sampling/60))
d_std = sqrt(diag(data{1}.Cd));
txt2 = sprintf('%s_2d_REJ%d_N%d_Nd%d_Nlu%d_Nr%d_%g_%g',txt,useRejection,N,Nd,Nlu,Nr,1000*d_std(1),1000*d_std(2));
%save(space2char(txt2,'_','\.'),'M_*','P_*', 'R_*')
%save(space2char(txt2,'_','\.'),'M_*','P_*', 'R_*','Mpr_*','Ppr_*', 'Rpr_*','D_*')
%save([space2char(txt2,'_','\.'),'_all'])

%%
figure(5);set_paper('landscape');clf;
subplot(1,3,1);imagesc(x,y,M_v_clay);title('v_{clay}');caxis([0 0.4]);axis image;colorbar
subplot(1,3,2);imagesc(x,y,M_sat_g);title('s_{gas}');axis image;caxis([0 1]);colorbar
subplot(1,3,3);imagesc(x,y,M_sat_o);title('s_{oil}');axis image;caxis([0 1]);colorbar
colormap(jet)
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'FontSize',8)
set(allAxesInFigure,'ydir','normal')
%axis(allAxesInFigure,'off')
for iax=1:length(allAxesInFigure);set(get(allAxesInFigure(iax),'xlabel'),'String','Inline','FontSize',10);end
for iax=1:length(allAxesInFigure);set(get(allAxesInFigure(iax),'ylabel'),'String','Crossline','FontSize',10);end
print_mul(sprintf('%s_post_mean_manus',txt2))


%%
Mpr_sat_g=reshape(mean(D_sat_g_log'),ny,nx);
Mpr_sat_o=reshape(mean(D_sat_o_log'),ny,nx);
Mpr_sat_b=reshape(mean(1-(D_sat_g_log'+D_sat_o_log')),ny,nx);;
Mpr_v_clay=reshape(mean(D_v_clay_log'),ny,nx);


figure(6);set_paper('landscape');clf;
subplot(1,3,1);imagesc(x,y,Mpr_v_clay);title('v_{clay}');caxis([0 0.4]);axis image;colorbar
subplot(1,3,2);imagesc(x,y,Mpr_sat_g);title('s_{gas}');axis image;caxis([0 1]);colorbar
subplot(1,3,3);imagesc(x,y,Mpr_sat_o);title('s_{oil}');axis image;caxis([0 1]);colorbar
colormap(jet)
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'FontSize',8)
set(allAxesInFigure,'ydir','normal')
%axis(allAxesInFigure,'off')
for iax=1:length(allAxesInFigure);set(get(allAxesInFigure(iax),'xlabel'),'String','Inline','FontSize',10);end
for iax=1:length(allAxesInFigure);set(get(allAxesInFigure(iax),'ylabel'),'String','Crossline','FontSize',10);end
print_mul(sprintf('%s_prior_mean_manus',txt2))

