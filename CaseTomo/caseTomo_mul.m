%% RAY
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='ray';is_slowness=1;caseTomo_setup

%%
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=500000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=500000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=500000;caseTomo_lookup_ml
%
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=500000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=500000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo1.mat';N=500000;caseTomo_lookup_ml


%% EIKONAL
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';is_slowness=1;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=0;caseTomo_setup
clear all;useCase='Kallerup';dx=0.1;forward.type='eikonal';is_slowness=1;caseTomo_setup

%%
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=500000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=5000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo0.mat';N=5000000;caseTomo_lookup_ml
%
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';di_use=1;caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=500000;di_use=1;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=5000000;di_use=1;caseTomo_rejection
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0_slo1.mat';N=5000000;caseTomo_lookup_ml


