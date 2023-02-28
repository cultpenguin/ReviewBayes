%% THE FINAL SETUP 

%% Setup / Parameterization
close all
% The linear forward
clear all;useCase='Kallerup';dx=0.10;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;caseTomo_setup
% The non-linear forward
clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;caseTomo_setup

%% The Inversions
% This first one is for Mina
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=200001;caseTomo_lookup_ml
% Sampling
% di_use: Use every di_use data
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=250000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';N=250000;di_use=1;caseTomo_rejection
% Linear Least Squares
clear all;fmat='caseTomo_Kallerup_dx10_Fray-ray_ME1_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-ray_ME1_slo0.mat';di_use=1;caseTomo_LeastSquares

return

%% OTHER TESTS

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
