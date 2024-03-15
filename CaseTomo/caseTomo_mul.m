
%% Setup / Parameterization
close all
% The linear forward
%clear all;useCase='Kallerup';dx=0.10;forward.type='ray';is_slowness=0;addStaticErr=1;caseTomo_setup
%clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;addStaticErr=1;caseTomo_setup
% The non-linear forward
clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';is_slowness=0;addStaticErr=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;addStaticErr=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';is_slowness=0;addStaticErr=2;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;addStaticErr=2;caseTomo_setup

clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';is_slowness=0;addStaticErr=3;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;addStaticErr=3;caseTomo_setup

clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';is_slowness=0;addStaticErr=3;doPriorGaussian=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;addStaticErr=3;doPriorGaussian=1;caseTomo_setup


clear all;useCase='Kallerup';dx=0.10;forward.type='fat';is_slowness=1;addStaticErr=3;caseTomo_setup
%% The inversions

clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE1_G0.mat';N=1000000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE2_G0.mat';N=1000000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE1_G0.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE2_G0.mat';N=1000000;di_use=1;caseTomo_rejection
%%

clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G0.mat';N=1000000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G0.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G1.mat';N=1000000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G1.mat';N=1000000;di_use=1;caseTomo_rejection

%% Localize rej for paper
close all
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0_SE3_G0.mat';N=10000;di_use=1;caseTomo_rejection_local

%% The Inversions
% This first one is for Mina
%clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=500001;caseTomo_lookup_ml
% 
% Arrays have incompatible sizes for this operation.
% 
% Error in sippi_abc_ml_predict (line 27)
%     D=(D-repmat(ABC.ml{iml}.mD',[1,size(D,2)]))./repmat(ABC.ml{iml}.vD',[1,size(D,2)]);
% 
% Error in caseTomo_lookup_ml (line 79)
% M_est=sippi_abc_ml_predict(ABC,D);
% 
% Error in caseTomo_mul (line 14)
% clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=500001;caseTomo_lookup_ml
% 
% Related documentation
% 

clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=4;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=4;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=4;caseTomo_rejection
return
% Sampling
% di_use: Use every di_use data
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=1000000;di_use=1;caseTomoc_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=1000000;di_use=4;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=1000000;di_use=4;caseTomo_rejection

% Linear Least Squares
clear all;fmat='caseTomo_Kallerup_dx10_Fray-ray_ME1_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-ray_ME1_slo0.mat';di_use=1;caseTomo_LeastSquares

return


%% TEST
clear all;close all
useCase='Kallerup';dx=0.10;forward.type='ray';is_slowness=0;
addStaticErr=2;
caseTomo_setup
%%
close all
fmat=txt;
%fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat'
N=100000;di_use=1;caseTomo_metropolis

%%
clear all;useCase='Kallerup';dx=0.10;forward.type='ray';is_slowness=0;caseTomo_setup
%%

%%
ccS=3;
ccR=3;
ccU=0;
[Ct, Cs, Cr, C]=correlated_traveltime_tomography_errors(forward.sources,forward.receivers,ccS,ccR,ccU);

Ct2= Ct+data{1}.Ct;

subplot(1,2,1)
imagesc(Ct2)
axis image
subplot(1,2,2)
d_err = gaussian_simulation_cholesky(0,Ct2+.0000001*eye(nd),1);plot(d_err)

%% OTHER TESTS

%% THE FINAL SETUP 

close all;clear all;
useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=0;doPriorGaussian=1;
caseTomo_setup
%clear all;
%fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';
%fmat=txt;
N=1000000;di_use=1;caseTomo_metropolis

return

% % RAY
% clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=0;caseTomo_setup
% clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';;di_use=1;caseTomo_LeastSquares
% clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';
% skipBadData=0;
% addExtraCt=0;
% N=100000;
% di_use=1;caseTomo_metropolis
% 
% % Eikonal
% clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=0;
% skipBadData=1;skipBadData=0;
% addExtraCt=0; % As we use the Eikonal solution!!!
% Do_comp_model_error=1-addExtraCt;
% caseTomo_setup
% %clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
% clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';N=100000-2;di_use=1;caseTomo_metropolis
% return
%% RAY
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=1;caseTomo_setup

%%
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';N=500000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';N=5000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';N=1000000;caseTomo_lookup_ml
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=5000000;caseTomo_lookup_ml
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
%
%clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
%clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=100000;di_use=1;caseTomo_metropolis
%clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=1000000;di_use=1;caseTomo_rejection
%clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=1000000;caseTomo_lookup_ml
%clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares


%% EIKONAL
clear all;useCase='Kallerup';dx=0.50;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=0;caseTomo_setup
%clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=1;caseTomo_setup
%clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=1;caseTomo_setup

%%
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=500000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=5000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=100000;caseTomo_lookup_ml
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-ray_ME1_slo0.mat';N=5000000;caseTomo_lookup_ml
clear all;fmat='caseTomo_Kallerup_dx50_Feikonal-ray_ME1_slo0.mat';N=5000000;caseTomo_lookup_ml                
%
%clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=100000;di_use=1;caseTomo_metropolis
%clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=1000000;di_use=1;caseTomo_rejection
%clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=1000000;caseTomo_lookup_ml


return

%%
for j=100:1:200;
C=data{1}.Ct(j,:);
cmin=min(C);
cmax=max(C);
nc=100;
cmap=jet(nc);
cmap = cmap_linear([1 1 1;1 1 .7;1 0 0;0 0 0],[0 .8 .9 1],nc);
colormap(cmap)
figure(99);clf
for i=1:size(data{1}.Ct,1);
    ic=ceil(interp1([cmin cmax],[1,nc],[C(i)]));
    x=[forward.sources(i,1),forward.receivers(i,1)];
    y=[forward.sources(i,2),forward.receivers(i,2)] ;
    lw=.01+2*(C(i)-cmin);
    plot(x , y,'k.')
    if ic>(nc*0.6)
        plot(x , y,'k','LineWidth',lw,'Color',cmap(ic,:));
    end
    hold on;
end
axis image
set(gca,'ydir','reverse')
drawnow;
end

%% TESTING
%%
clear all;close all;
useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;
doDataGaussian=0;doDataBoost=1;useCorrelatedNoise=0;skipBadData=0;
caseTomo_setup
fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=50000;di_use=1;
caseTomo_metropolis

clear all;close all;
useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;
doDataGaussian=0;doDataBoost=1;useCorrelatedNoise=1;skipBadData=0;
caseTomo_setup
fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=50001;di_use=1;
caseTomo_metropolis

clear all;close all;
useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;
doDataGaussian=0;doDataBoost=0;useCorrelatedNoise=0;skipBadData=0;
caseTomo_setup
fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=50002;di_use=1;
caseTomo_metropolis


clear all;close all;
useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;
doDataGaussian=0;doDataBoost=0;useCorrelatedNoise=1;skipBadData=0;
caseTomo_setup
fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=50003;di_use=1;
caseTomo_metropolis

% Gaussian
clear all;close all;
useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;
doDataGaussian=1;doDataBoost=1;useCorrelatedNoise=1;skipBadData=0;
caseTomo_setup
fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=50010;di_use=1;
caseTomo_metropolis

clear all;close all;
useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;
doDataGaussian=1;doDataBoost=0;useCorrelatedNoise=1;skipBadData=0;
caseTomo_setup
fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=50011;di_use=1;
caseTomo_metropolis

caseTomo_LeastSquares
