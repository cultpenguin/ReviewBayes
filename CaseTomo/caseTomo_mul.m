%% Setup 
rng('default')
%clear all;dx=0.1;forward.type='eikonal';caseTomo_setup
%clear all;dx=0.25;forward.type='eikonal';caseTomo_setup
%clear all;dx=0.5;forward.type='eikonal';caseTomo_setup
%clear all;dx=0.1;forward.type='fat';caseTomo_setup
%clear all;dx=0.25;forward.type='fat';caseTomo_setup
%clear all;dx=0.5;forward.type='fat';caseTomo_setup

%% METROPOLIS
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=200000;caseTomo_metropolis

clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=200001;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=200001;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=200001;caseTomo_metropolis

%% Rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection

