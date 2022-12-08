%% Setup 

clear all;dx=0.1;forward.type='fat';caseTomo_setup
clear all;dx=0.25;forward.type='fat';caseTomo_setup
clear all;dx=0.5;forward.type='fat';caseTomo_setup

%%
%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=50001;caseTomo_rejection
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=50001;caseTomo_rejection
%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=50000;caseTomo_rejection
%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=50001;caseTomo_rejection
%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=50002;caseTomo_rejection

%clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=20000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=20000;caseTomo_metropolis
%clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=20000;caseTomo_metropolis
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=40001;caseTomo_metropolis
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=40001;caseTomo_metropolis
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=40001;caseTomo_metropolis

%% Rejection

%fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=1000;caseTomo_rejection
%fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=100000;caseTomo_rejection
%fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=1000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000000;caseTomo_rejection
clear all;fmat='caseBayesian_dx50_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx25_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection
clear all;fmat='caseBayesian_dx10_Ffat-none_ME0.mat';N=5000001;caseTomo_rejection

