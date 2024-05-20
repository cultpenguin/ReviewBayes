%%
delete('case*mat')

progress_out('Start case 2, 2D',1)
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

