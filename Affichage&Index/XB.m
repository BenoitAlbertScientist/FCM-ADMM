function [xb,compacity,separability]=XB(X,U,V,parameters)
% Xie-Beni XB and XBMW
%    [xb,compacity,separability]=XB(x,U,V,m,parameters)
% 
% INPUTS
%   X: input matrix (objects x attributs)
%   U: partition matrix  (objects x cluster)
%   V: centroids    (cluster x attributs)
%   Optional :  
%       parameters.choice_index : 0 =  XB (default) 
%                                 1 =  XBMW
%
%       parameters.give_cov      : #For Mahalanobis
%                                  0 = give cov matrix (default), 1=give inverse cov
%       parameters.matrix        : martrix of clusters 
%                                    (cluster x (attributs x attributs))
%       parameters.iprint        : 0 = no print (default)
%       parameters.m             : fuzzy parameter = 2 (default) 
%
% OUTPUTS
%   xb           : Xie Beni index
%   compacity    : up part 
%   separability : low part
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 06-21-2023
%  version: 1
%  --------------------------------------------------------------------------

% ------------------------ Check parameters -------------------------
if ~isfield(parameters,'choice_index') parameters.choice_index=0;end
if ~isfield(parameters,'give_cov') parameters.give_cov=0;end
if ~isfield(parameters,'iprint') parameters.iprint=0;end
if ~isfield(parameters,'m') parameters.m=2;end

if (parameters.choice_index~=0 && parameters.choice_index~=1) parameters.choice_index=0;end
if (parameters.give_cov~=0 && parameters.give_cov~=1) parameters.give_cov=0;end
if (parameters.iprint~=0 && parameters.iprint~=1) parameters.iprint=0;end


choice_index=parameters.choice_index;
give_cov=parameters.give_cov;
iprint=parameters.iprint;
m=parameters.m;

[n,c] = size(U);[cv,nd] = size(V);[nx,ndx] = size(X);

if n~=nx
    disp('Wrong number of objects : U (objects x clusters) X(objects x attributs)')
end

if c~=cv
    disp('Wrong number of clusters : U (objects x clusters) V(clusters x attributs)')
end

if nd~=ndx
    disp('Wrong number of attributs : X (objects x attributs) V(clusters x attributs)')
    
end

% ------------------------ Init -------------------------


if iprint == 1
    nom = 'XB'; cov_or_inv={'','INVERSE of'};
    if choice_index==0
        nom = strcat(' XBMW with',cov_or_inv{give_cov},'fuzzy covariance : \n');
    end
    fprintf('*******************************************\n');
    disp(nom);
    fprintf('-------------------------------------------\n');
end

compacity = 0;
separability = 0;
Cov=cell(c,1);S=cell(c,1);

%Cov et InvCov
if(choice_index==0)
    for j=1:c
        Cov{j}=eye(nd);
        S{j}=eye(nd);
    end
else
    if (give_cov==0)
        Cov=parameters.matrix;
        for j=1:c
            S{j}=inv(Cov{j});
        end
    else
        S=parameters.matrix;
        for j=1:c
            Cov{j}=inv(S{j});
        end 
    end
end

%Compacity
for j=1:c
    for i=1:n
        dij=(X(i,:)-V(j,:))*S{j}*(X(i,:)-V(j,:))';
        compacity = compacity +U(i,j)^m*dij;
    end
end

%Separbility
dv=Inf(c);
for j1=1:(c-1)
    for j2=(j1+1):c
         dv(j1,j2)=norm(V(j1,:)-V(j2,:))^2+...
             trace(Cov{j1}+Cov{j2}-2*sqrtm(sqrtm(Cov{j1})*Cov{j2}*sqrtm(Cov{j1})));
    end
end

separability=min(min(dv));
xb=compacity/(n*separability);

if iprint == 1
    fprintf(num2str(xb,'%1.2e'));
    fprintf('-------------------------------------------\n');
end

end
