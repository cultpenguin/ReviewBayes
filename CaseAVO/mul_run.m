%%
%clear all;close all;useCase=2;Nlu=50000;Nlsr=100;CaseAVO_inv_ABC_2d
%clear all;close all;useCase=3;Nlu=50000;Nr=100;CaseAVO_inv_ABC_2d
%clear all;close all;useCase=4;Nlu=50000;Nr=100;CaseAVO_inv_ABC_2d
delete('case*mat')
progress_out('Start case 3, eikonal',1)
clear all;close all;useCase=3;Nlu=50000;Nr=400;useRejection=1;CaseAVO_inv_ABC
progress_out('Start case 4, eikonal')
clear all;close all;useCase=4;Nlu=50000;Nr=400;useRejection=1;CaseAVO_inv_ABC
progress_out('Start case 2, eikonal')
clear all;close all;useCase=2;Nlu=50000;Nr=400;useRejection=1;CaseAVO_inv_ABC
progress_out('Start case 3, eikonal')

