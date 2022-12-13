%% Setup 
rng('default')
%clear all;dx=0.1;forward.type='eikonal';caseTomo_setup
%clear all;dx=0.25;forward.type='eikonal';caseTomo_setup
%clear all;dx=0.5;forward.type='eikonal';caseTomo_setup
%clear all;dx=0.1;forward.type='fat';caseTomo_setup
%clear all;dx=0.25;forward.type='fat';caseTomo_setup
%clear all;dx=0.5;forward.type='fat';caseTomo_setup
%clear all;dx=0.1;forward.type='ray';caseTomo_setup
%clear all;dx=0.25;forward.type='ray';caseTomo_setup
%clear all;dx=0.5;forward.type='ray';caseTomo_setup

%% Least Squares
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';caseTomo_LeastSquares

clear all;fmat='caseBayesian_dx50_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx25_Fray-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx10_Fray-none_ME0.mat';caseTomo_LeastSquares



%%  small test
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=20000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200000;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';caseTomo_LeastSquares
clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=20000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=200000;caseTomo_rejection
%clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';caseTomo_LeastSquares

%% METROPOLIS
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx50_Feikonal-none_ME0.mat';N=200000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx25_Feikonal-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Feikonal-none_ME0.mat';N=200000;caseTomo_metropolis

%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=200001;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200001;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=200001;caseTomo_metropolis

%% Rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection

