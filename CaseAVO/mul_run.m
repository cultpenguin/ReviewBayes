progress_out('Start CaseAVO',1)

%%
progress_out('Start case wx0,dx1 localprior')
clear all;close all;wx=0;dx=1;mina_load_localprior_N;
progress_out('Start case wx1,dx1 localprior')
clear all;close all;wx=1;dx=1;mina_load_localprior_N;
progress_out('Start case wx2,dx1 localprior')
clear all;close all;wx=2;dx=1;mina_load_localprior_N;
progress_out('Start case wx3,dx1 localprior')
clear all;close all;wx=3;dx=1;mina_load_localprior_N;


return

%%
delete('case*mat')

progress_out('Start case 2, 2D')
clear all;close all;useCase=2;Nlu=50000;Nr=100;CaseAVO_inv_ABC_2d
progress_out('Start case 3, 2D')
clear all;close all;useCase=3;Nlu=50000;Nr=100;CaseAVO_inv_ABC_2d
progress_out('Start case 4, 2D')
clear all;close all;useCase=4;Nlu=50000;Nr=100;CaseAVO_inv_ABC_2d

progress_out('Start case 3')
clear all;close all;useCase=3;Nlu=50000;Nr=400;useRejection=1;CaseAVO_inv_ABC
progress_out('Start case 4')
clear all;close all;useCase=4;Nlu=50000;Nr=400;useRejection=1;CaseAVO_inv_ABC
progress_out('Start case 2')
clear all;close all;useCase=2;Nlu=50000;Nr=400;useRejection=1;CaseAVO_inv_ABC
progress_out('Start case 3, eikonal')

%%
delete('case*mat')
delete('MINA*mat')
progress_out('START mina_load')
mina_load
progress_out('END mina_load')



