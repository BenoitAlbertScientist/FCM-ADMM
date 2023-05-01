%Téléchargement des données
addpath 'Affichage&Index'
addpath 'Data'
load iris_n

n=size(x,2); nd=size(x,1);
c=length(cl);


%% APPLICATION DE FCM
rng('default'); %INITIALISATION DE L'aléatoire


parameters.init = 1;
parameters.distance = 1;
parameters.iprint = 1;

%ADMM
parameters.ncadmm = 5;
parameters.r = 30;
[u,v,S,iter,fobj] = FCM_ADMM(x,c,parameters);
fprintf("ARI  = %10.4f \n",ARI(HP,Fuzzy2Hard(u)));
titre = strcat(name_data,'[ADMM]');
DisplayClustering2D(x',v',Fuzzy2Hard(u),S,titre);

%AO
[u,v,S,iter,fobj] = FCM_AO(x,c,parameters);
fprintf("ARI  = %10.4f \n",ARI(HP,Fuzzy2Hard(u)));
titre = strcat(name_data,'[AO]');
DisplayClustering2D(x',v',Fuzzy2Hard(u),S,titre);   
  
