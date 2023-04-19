%% Setup ABC
if ~exist('txt','var')
    txt = 'MINA_ABC';
end

mat_file = sprintf('%s_N%d',txt,Nlu);

%% Setup ABC structure
disp('Setup ABC')
clear ABC
Nlu=min([N Nlu]);
i_use=randomsample(1:N,Nlu);

if exist([mat_file,'.mat'])==2
    disp(sprintf('Loading %s',mat_file))
    load(mat_file)
else

    i=0;
    for k=i_use;1:1000000;N;5000;%N;
        i=i+1;
        if mod(i,1000)==0;progress_txt(i,length(i_use));end
        % m
        useMstyle=1;
        if useMstyle==1;
            j=0;
            j=j+1;ABC.m{i}{j}=D.sat_g_log(k);
            j=j+1;ABC.m{i}{j}=D.sat_o_log(k);
            j=j+1;ABC.m{i}{j}=D.v_clay_log(k);
            j=j+1;ABC.m{i}{j}=D.v_p(k);
            j=j+1;ABC.m{i}{j}=D.v_s(k);
            j=j+1;ABC.m{i}{j}=D.rho(k);
            j=j+1;ABC.m{i}{j}=D.depth_all(k);
            %j=j+1;ABC.m{i}{j}=D.porosity(k);
        else
            ABC.m{i}{1}=[D.sat_g_log(k) D.sat_o_log(k) D.v_clay_log(k) D.v_p(k) D.v_s(k) D.rho(k) D.depth_all(k)];
        end

        if i==1;
            ABC.m=repmat(ABC.m,[length(i_use),1]);
        end

        % d
        j=0;
        %j=j+1;ABC.d{i}{j}=[D.g_sim(k), D.r_0_sim(k)];
        %j=j+1;ABC.d{i}{j}=[D.depth_all(k)];
        j=j+1;ABC.d{i}{j}=[ D.r0_sim_no_noise(k), D.g_sim_no_noise(k), D.depth_all(k)]';
        j=j+1;ABC.d{i}{j}=[ D.r_0_sim(k), D.g_sim(k), D.depth_all(k)]';
        
        if i==1;
            ABC.d=repmat(ABC.d,[length(i_use),1]);
        end

    end
    disp(sprintf('Saving %s',mat_file))
    save(mat_file,'ABC');
end

Nlu=length(ABC.m);