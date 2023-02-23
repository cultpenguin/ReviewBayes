

%% RAY
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=1;caseTomo_setup

%%
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=100000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=1000000;caseTomo_lookup_ml
%
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=100000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=1000000;caseTomo_lookup_ml


%% EIKONAL
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=1;caseTomo_setup

%%
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=100000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=1000000;caseTomo_lookup_ml
%
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=100000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=1000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=1000000;caseTomo_lookup_ml


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
