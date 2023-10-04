function [u,v,S,iter,fobj] = FCM_AO(x,c,parameters)
% Fuzzy C-Means with Alternating Optimization
%    [u,v,S,iter,fobj] = FCM_AO(x,c,parameters)
% 
% INPUTS
%   x: input matrix p*n (attributs x objects)
%   c: number of desired clusters
%   Optional : 
%       parameters.distance: 0=euclidean, 1=mahalanobis
%       parameters.init    : 0=random initialization of the center, 
%                            1=initialization with
%                              ADMMeuclidian distance=0,itmax=50,r=2.5.
%                            2=specific initialization (gives centers g)
%       parameters.g
%       parameters.tol     : tolerance (default:10^-3)
%       parameters.itmax   : maximal number iterations maximal (default:1000)
%       parameters.
%           iprint         : 1=display general informations (default:1,else:0)
%           iprint_inside  : 1=display at each iteration (default:0)
%
%
% OUTPUTS
%   u: fuzzy partition  (clusters x objects)
%   v: centroids (attributs x clusters)
%   S: cell of cov-matrix (clusters (attributs x attributs ))
%   iter: number of iterations
%   fobj : objectif function
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 10-01-2023
%  version: 2
%  --------------------------------------------------------------------------


% imensions
if nargin<2  
    error('FCM needs two arguments');
else
    n=size(x,2); nd=size(x,1);
end


% ------------------------ Check parameters -------------------------
if ~isfield(parameters,'init') parameters.init=1;end
if ~isfield(parameters,'distance') parameters.distance=0;end
if ~isfield(parameters,'tol') parameters.tol=10^-3;end
if ~isfield(parameters,'itmax') parameters.itmax=1000;end
if ~isfield(parameters,'iprint') parameters.iprint=1;end
if ~isfield(parameters,'iprint_inside') parameters.iprint_inside=0;end

if (parameters.init~=0 && parameters.init~=1 && parameters.init~=2) parameters.init=1;end
if (parameters.distance~=0 && parameters.distance~=1) parameters.distance=0;end
if parameters.itmax<1 parameters.itmax=1000;end
if (parameters.iprint~=0 && parameters.iprint~=1) parameters.iprint=1;end
if (parameters.iprint_inside~=0 && parameters.iprint_inside~=1) parameters.iprint_inside=0;end

init=parameters.init;
distance=parameters.distance;
tol=parameters.tol;
itmax=parameters.itmax;
iprint=parameters.iprint;
iprint_inside=parameters.iprint_inside;

if iprint == 1
    fprintf('*******************************************\n');
    fprintf('\t Fuzzy C-means with AO\n');
    fprintf('-------------------------------------------\n');
    fprintf('Number of objects  = %5i\n',n);
    fprintf('Number of clusters = %5i\n',c);
end

% ---------------------- Initialization ----------------------

eps0=10^-10; %avoid fuzzy covariance to be null 


if init == 0 
% Random
    % Mass
    u=rand(n,c); su=sum(u,2);
    u=u./su(:,ones(1,c));

    % Center of mass
    u2=u.^2; su=sum(u2); v=x*u2;
    v=v./su(ones(nd,1),:);
    if iprint == 1;fprintf('Initialization : Random\n');end;
    
elseif init == 1
% Initialization with Euclidean distance    
    parameters_euc.init = 0;
    parameters_euc.distance = 0;
    parameters_euc.r = 2.5;
    parameters_euc.itmax = 50;
    parameters_euc.iprint = 0;
    [u_eu,v_eu,S_eu,iter_eu] = FCM_ADMM(x,c,parameters_euc);
    %Remarque : Initialization could be FCM_AO with euclidian
    u = u_eu; v = v_eu; S = S_eu;
    if iprint == 1;fprintf('Initialization : Euclidan[iter=%2i  (max. 50)]\n',iter_eu);end;

elseif init == 2
% Specific initialization
    v = parameters.g;
    for i=1:n
        su=0;
        for k=1:c
            dd=(x(:,i)-v(:,k))'*(x(:,i)-v(:,k));
            su=su+1/dd;
        end
        for j=1:c
            dd=(x(:,i)-v(:,j))'*(x(:,i)-v(:,j));
            u(i,j)=max(0,1/dd/su);
        end
    end
    if iprint == 1;fprintf('Initialization : Specific\n');end;
    
end

% Initialization of Var-Covariance
eps0 = 10^-8;
if distance == 0
    for j=1:c;S{j}=eye(nd);end 
else 
    ux=u(:)';  ic=[1:c]'; ic=kron(eye(c),ones(1,n))'*ic;
    d=repmat(x,1,c)-v(:,ic');
    p=ux(ones(nd,1),:).*d;
    for j=1:c
        Sj=eps0*eye(nd); j1=(j-1)*n+1; j2=(j-1)*n+n;
        S{j}=p(:,j1:j2)*p(:,j1:j2)'+eps0*eye(nd);
        scal=(det(S{j}))^(1/nd); 
        S{j}=scal*inv(Sj);
    end 
end


% ------------------------- Iterations ----------------------------------
err=1;
iter=0;

while (iter < itmax && err>tol)
    iter=iter+1;
    uu=u; vv=v; SS=S;
    
    % Compute centers : V
    su=sum(u.*u,1);
    for j=1:c
        v(:,j)=zeros(nd,1);
        for i=1:n
            v(:,j)=v(:,j)+u(i,j)*u(i,j)*x(:,i);
        end
        v(:,j)=v(:,j)/su(j);
    end
    
    % Compute distance matrix : S
    if distance == 1
        for j=1:c
            Sj=eps0*eye(nd);

            su=0;
            for i=1:n
                Sj=Sj+u(i,j)*u(i,j)*(x(:,i)-v(:,j))*(x(:,i)-v(:,j))';
            end
            dS(j)=det(Sj); 
            scal=dS(j)^(1/nd);
            S{j}=scal*inv(Sj);
        end 
    end

    % Compute partition :  U
    for i=1:n
        su=0;
        for k=1:c
            dd=(x(:,i)-v(:,k))'*S{k}*(x(:,i)-v(:,k))+eps0;
            su=su+1/dd;
        end
        for j=1:c
            dd=(x(:,i)-v(:,j))'*S{j}*(x(:,i)-v(:,j))+eps0;
            u(i,j)=max(0,1/dd/su);
        end
    end
    
    
    % Compute errors (stopping criterion)
    ndu=sum(sum((u-uu).^2));   nu=sum(sum(u.^2));
    ndv=sum(sum((v-vv).^2));   nv=sum(sum(v.^2));
     
    ndS=0; nS=0;
    for j=1:c
        nS=nS+norm(S{j},'fro')^2;
        ndS=ndS+norm(S{j}-SS{j},'fro')^2;
 
    end
    errd=[ndu;ndv;ndS];
    errn=[nu;nv;nS];
    
    err=sqrt(sum(errd)/sum(errn));

    if (iprint_inside == 1)    
        %Objectif function
        fobj = 0;
        for j=1:c
           Sj = S{j};
           for i=1:n
              pij = u(i,j)*(x(:,i)-v(:,j));
              fobj = fobj + pij'*Sj*pij;
           end
        end
       fprintf(">>iter=%i | err=%1.6e\n | J_FCM=%15.8e\n",iter,err,fobj);
    end
 
    
end

%Function
fobj = 0;
for j=1:c
   Sj = S{j};
   for i=1:n
      pij = u(i,j)*(x(:,i)-v(:,j));
      fobj = fobj + pij'*Sj*pij;
   end
end

if (iprint == 1)   
    fprintf('-------------------------------------------\n');
    fprintf("Objectif function F(U,V,S) =%e\n",fobj);
    fprintf("[iter=%i (max. %i)| err=%1.6e]\n",iter,itmax,err);
    fprintf('*******************************************\n');
   
end

end
