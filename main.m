%% Download data 
addpath 'Affichage&Index'
addpath 'Data'
load iris_n

n=size(x,2);    %Number of objects
nd=size(x,1);   %Number of attributs
c=length(cl);   %Number of clusters

%% COMPARAISON ADMM vs AO
%Apply on FCM-GK model. 

parameters.init = 1; %Init 
parameters.distance = 1; %Mahalanobis distance
parameters.iprint = 1;
%ADMM
name_meth = 'ADMM'; rng('default'); %Rand init
parameters.ncadmm = 5;parameters.r = choix_r_ADMM(name_data);
parameters.tol = 10^-4;
[u,v,S,iter,fobj] = FCM_ADMM(x,c,parameters);
EVAL(x,u,v,S,HP,name_data,name_meth);

% ----- AO
name_meth = 'AO'; rng('default'); %Rand init
parameters.tol = 10^-3;
[u,v,S,iter,fobj] = FCM_AO(x,c,parameters);
EVAL(x,u,v,S,HP,name_data,name_meth);

%% Evaluation 
%Evaluation with ARI, PE, XB and XBMW.
%Print in 2D clustering.

function [] = EVAL(x,u,v,S,HP,name_data,name_meth)

    %ARI
    hp=Fuzzy2Hard(u);
    fprintf("ARI  = %.2f \n",ARI(HP,Fuzzy2Hard(u)));
    
    %PE
    fprintf("PE  = %1.2f \n",PE(u));
    
    %FS
    fprintf("FS  = %1.2f \n",FS(x',u));
    
    %XB
    parameters_XB.choice_index=0;
    xb=XB(x',u,v',parameters_XB);
    fprintf("XB  = %1.2f \n",xb);
    
    %XBMW
    parameters_XBMW.choice_index=1;
    parameters_XBMW.give_cov=1; %S is inverse of covariance matrix
    parameters_XBMW.matrix=S;
    xbmw=XB(x',u,v',parameters_XBMW);
    fprintf("XBMW  = %1.2f \n",xbmw);
    
    %DISPLAY
    titre = strcat(name_data,name_meth);
    DisplayClustering2D(x',v',hp,S,titre);

end

%% FIX r ADMM

function [r] = choix_r_ADMM(name)
    if strcmp(name,'iris_n');r=320;
    else;if strcmp(name,'algfor_n');r=10^6; 
    else;if strcmp(name,'drybean_n');r=10^5; 
    else;if strcmp(name,'glass2_n');r=10^5;
    else;if strcmp(name,'seed_n');r=10^5;
    else;if strcmp(name,'letterIJL_n');r=10^4;
    else;if strcmp(name,'wdbc_n');r=10^5; 
    else;if strcmp(name,'wifi_n');r=250;
    else;if strcmp(name,'wine_n');r=400;
   
    else;if strcmp(name,'a1_n');r=40;
    else;if strcmp(name,'a3_n');r=40;
    else;if strcmp(name,'asymmetric_n');r=10;
    else;if strcmp(name,'dim032_n');r=1000;
    else;if strcmp(name,'dim064_n');r=500;
    else;if strcmp(name,'s1_n');r=10^3;
    else;if strcmp(name,'s3_n');r=10^3;
    else;if strcmp(name,'skewed_n');r=10;
    else; disp(strcmp('choix_r_ADMM ',name,' unknows'));
    end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end;end

end
