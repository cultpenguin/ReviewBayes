%clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';caseTomo_setup
%clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat';N=200001;ml.normalize=1;caseTomo_lookup
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat';N=20001;ml.normalize=1;caseTomo_lookup_ml
return
%%
clear all;useCase='Kallerup';dx=0.25;forward.type='ray';caseTomo_setup
clear all;useCase='Kallerup';dx=0.25;forward.type='eikonal';caseTomo_setup
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat';N=250001;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0.mat';N=250001;caseTomo_metropolis

clear all;useCase='Kallerup';dx=0.10;forward.type='ray';caseTomo_setup
clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';caseTomo_setup
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0.mat';N=250001;caseTomo_metropolis
clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0.mat';N=250001;caseTomo_metropolis
return

%% Setup 

clear all;useCase='Kallerup';dx=0.05;forward.type='ray';caseTomo_setup
clear all;useCase='Kallerup';dx=0.10;forward.type='ray';caseTomo_setup
clear all;useCase='Kallerup';dx=0.10;forward.type='fat';caseTomo_setup
clear all;useCase='Kallerup';dx=0.10;forward.type='eikonal';caseTomo_setup
clear all;useCase='Arrenaes';dx=0.10;forward.type='ray';caseTomo_setup
clear all;useCase='Arrenaes';dx=0.10;forward.type='fat';caseTomo_setup
clear all;useCase='Arrenaes';dx=0.10;forward.type='eikonal';caseTomo_setup

%%
%clear all;fmat='caseTomo_Kallerup_dx5_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0.mat';N=50001;caseTomo_metropolis
%clear all;fmat='caseTomo_Kallerup_dx10_Ffat-none_ME0.mat';caseTomo_LeastSquares
%clear all;fmat='caseTomo_Kallerup_dx10_Ffat-none_ME0.mat';N=50000;caseTomo_metropolis
%clear all;fmat='caseTomo_Kallerup_dx10_Feikonal-none_ME0.mat';N=50000;caseTomo_metropolis
%%
return

return
%%
rng('default')
clear all;dx=0.10;forward.type='eikonal';caseTomo_setup
clear all;dx=0.25;forward.type='eikonal';caseTomo_setup
clear all;dx=0.50;forward.type='eikonal';caseTomo_setup
clear all;dx=0.10;forward.type='fat';caseTomo_setup
clear all;dx=0.25;forward.type='fat';caseTomo_setup
clear all;dx=0.50;forward.type='fat';caseTomo_setup
clear all;dx=0.10;forward.type='ray';caseTomo_setup
clear all;dx=0.25;forward.type='ray';caseTomo_setup
clear all;dx=0.50;forward.type='ray';caseTomo_setup

%%  small test
clear all;
useCase='Kallerup';
useCase='Arrenaes';
fmat='caseTomo_Kallerup_dx10_Ffat-none_ME0.mat';caseTomo_LeastSquares

%clear all;fmat='caseTomo_Kallerup_dx10_Fray-none_ME0.mat';caseTomo_LeastSquares
%clear all;fmat='caseTomo_Kallerup_dx25_Ffat-none_ME0.mat';N=20000;caseTomo_metropolis
%clear all;fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat';N=20000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=20000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200000;caseTomo_rejection
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';caseTomo_LeastSquares
%clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=20000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=200000;caseTomo_rejection
%clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';caseTomo_LeastSquares

%% Least Squares
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';caseTomo_LeastSquares

clear all;fmat='caseBayesian_dx50_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx25_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx10_Fray-none_ME0.mat';caseTomo_LeastSquares

%% METROPOLIS
clear all;fmat='caseBayesian_dx50_Fray-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Fray-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Fray-none_ME0.mat';N=200000;caseTomo_metropolis

clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis

clear all;fmat='caseBayesian_dx50_Feikonal-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Feikonal-none_ME0.mat';N=200000;caseTomo_metropolis

%% Rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection

clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection

clear all;fmat='caseBayesian_dx50_Feikonal-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Feikonal-none_ME0.mat';N=5000000;caseTomo_rejection
