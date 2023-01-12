# ReviewBayes

## Matlab packages

    git clone git@github.com:cultpenguin/sippi.git
    git clone git@github.com:cultpenguin/mgstat.git
    git clone git@github.com:cultpenguin/sippi-abc.git
  
## Set paths in Matlab

    addpath sippi
    addpath sippi-abc
    sippi_set_path
    
    


## Setting up some data

    clear all;
    useCase='Kallerup';
    dx=0.25;
    forward.type='eikonal';
    caseTomo_setup

This creates 'caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat' that contains the data, prior, forward amd noise model

Themn to create a large set of model and data (and data with noise) use 

    clear all;
    fmat='caseTomo_Kallerup_dx25_Feikonal-none_ME0.mat';
    N=200001;
    ml.normalize=1;
    
The simulated models and data are available in 
   
    ABC.m{i}{1} : Velocity model 1
    ABC.d{i}{1} : traveltime for model 1
    ABC.d{i}{2} : traveltime WITH noise for model 1
    
   
    
    

  
