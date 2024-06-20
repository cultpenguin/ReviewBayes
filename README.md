# ReviewBayes

## CaseAVO


## CaseTomo -- KALLERUP GPR DATA
The data are available through 
B. B. Jensen, T. M. Hansen, L. Nielsen, K. S. Cordua, N. Tuxen, A. Tsitonaki, M. C. Looms, <b>Accounting for Modeling Errors in Linear Inversion of Crosshole Ground-Penetrating Radar Amplitude Data: Detecting Sand in Clayey Till</b>.

[https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022JB024666]

### Matlab packages

    git clone git@github.com:cultpenguin/sippi.git
    git clone git@github.com:cultpenguin/mgstat.git
    git clone git@github.com:cultpenguin/sippi-abc.git
  
### Set paths in Matlab

    addpath sippi
    addpath mgstat
    addpath sippi-abc
    sippi_set_path
    
    


### Setting up some data

    clear all;
    useCase='Kallerup';
    dx=0.25;
    forward.type='eikonal';
    caseTomo_setup

This creates 'caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat' that contains the data, prior, forward amd noise model

Themn to create a large set of model and data (and data with noise) use 

    clear all;
    fmat='caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0.mat';
    N=200001;
    ml.normalize=1;
    caseTomo_lookup_ml       
    
The simulated models and data are available in the HDF5 file
caseTomo_Kallerup_dx10_Feikonal-ray_ME1_slo0_lu_noise_N200001.h5, as

    (base) tmeha@DellKontor:/mnt/f/PROGRAMMING/ReviewBayes/CaseTomo$ h5ls caseTomo_Kallerup_dx25_Feikonal-none_ME0_lu_N200001.h5
    D1                       Dataset {200001, 412} -- DATA WITHOUT NOISE
    D2                       Dataset {200001, 412} -- DATA WITH NOISE
    M1                       Dataset {200001, 551} -- VELOCITY MODEL
    M2                       Dataset {200001, 1}
    M3                       Dataset {200001, 551}
    (base) tmeha@DellKontor:/mnt/f/PROGRAMMING/ReviewBayes/CaseTomo$      
    ABC.m{i}{1} : Velocity model 1
    ABC.d{i}{1} : traveltime for model 1
    ABC.d{i}{2} : traveltime WITH noise for model 1
    
    x=prior{1}.x
    y=prior{1}.y


    
   
    
    

  
