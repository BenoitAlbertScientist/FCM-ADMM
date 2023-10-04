function [u,v,S,iter,fobj] = FCM_ADMM(x,c,parameters)
% Fuzzy C-Means with ADMM method
%    [u,v,S,iter] = FCM_ADMM(x,c,parameters)
% 
% INPUTS
%   x: input matrix p*n (attributs x objects)
%   c: number of desired clusters
%   Optional : 
%       parameters.distance: 0=euclidean, 1=mahalanobis
%       parameters.init    : 0=random initialization of the center, 
%                            1=initialization with
%                              ADMMeuclidian distance=0,itmax=50,r=2.5.
%                            2=specific initialization (gives centers)
%       parameters.g
%       parameters.r       : penalty parameter (default:4*attributs*objects*clusters)
%       parameters.ncadmm  : number of sub ADMM boucle (default:5 to insure convergence)
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
%   iter: numbers of iteration
%   fobj : objectif function
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 03-25-2023
%  version: 1
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
if ~isfield(parameters,'r') parameters.r=4*n*nd*c;end
if ~isfield(parameters,'ncadmm') parameters.ncadmm=5;end
if ~isfield(parameters,'tol') parameters.tol=10^-3;end
if ~isfield(parameters,'itmax') parameters.itmax=1000;end
if ~isfield(parameters,'iprint') parameters.iprint=1;end
if ~isfield(parameters,'iprint_inside') parameters.iprint_inside=0;end

if (parameters.init~=0 && parameters.init~=1 && parameters.init~=2) parameters.init=1;end
if (parameters.distance~=0 && parameters.distance~=1) parameters.distance=0;end
if parameters.r<1 parameters.r=4*n*nd*c;end
if parameters.ncadmm<1 parameters.ncadmm=5;end
if parameters.itmax<1 parameters.itmax=1000;end
if (parameters.iprint~=0 && parameters.iprint~=1) parameters.iprint=1;end
if (parameters.iprint_inside~=0 && parameters.iprint_inside~=1) parameters.iprint_inside=0;end

init=parameters.init;
distance=parameters.distance;
r=parameters.r;
ncadmm=parameters.ncadmm;%nb sub boucle in ADMM
tol=parameters.tol;
itmax=parameters.itmax;
iprint=parameters.iprint;
iprint_inside=parameters.iprint_inside;

if ~isfield(parameters,'HP') parameters.HP=0;end
HP=parameters.HP;


if iprint == 1
    fprintf('*******************************************\n');
    fprintf('\t Fuzzy C-means with ADMM\n');
    fprintf('-------------------------------------------\n');
    fprintf('Number of objects  = %5i\n',n);
    fprintf('Number of clusters = %5i\n',c);
    fprintf('r = %5i | sub ADMM boucle = %2i\n',r,ncadmm);
end

% ---------------------- Initialization ----------------------

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
    u = u_eu; v = v_eu; S = S_eu;
    if iprint == 1;fprintf('Initialization : Euclidan[iter=%2i  (max. 50)]\n',iter_eu);end;
elseif init == 2
%Specific
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

% Auxiliares variables and multipliers
ux=u(:)';  ic=[1:c]'; ic=kron(eye(c),ones(1,n))'*ic;
d=repmat(x,1,c)-v(:,ic');
p=ux(ones(nd,1),:).*d;
lam2=2*p;
lam1=ux(ones(nd,1),:).*lam2;

% Inv- Var-Covariance
eps0 = 10^-8;
if distance == 0
    for j=1:c
        S{j}=eye(nd);
    end 
else 
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
    pp=p; dd=d;
    
    % Inner iter to ensure convergence
    for k=1:ncadmm
        % Compute partition :  U
        lamp=lam2+r*p;
        ddij=sum(d.^2,1); dmij=sum(d.*lamp,1);
        ddij=reshape(ddij,n,c); dmij=reshape(dmij,n,c);
        ai=sum(1./(ddij+eps0),2); ai=ai(:,ones(1,c));
        bi=sum(dmij./(ddij+eps0),2); bi=bi(:,ones(1,c));
        u=(dmij-bi./ai+r./ai)./(ddij+eps0)/r;
        u=max(0,u);
 
        % Compute centers : V
        dl=d+lam1/r;
        v=zeros(nd,c);
        for j=1:c
            j1=(j-1)*n+1; j2=(j-1)*n+n;
            v(:,j)=sum(x-dl(:,j1:j2),2)/n;
        end
        
        % Compute distance matrix : S
        if distance == 1
            for j=1:c
                j1=(j-1)*n+1; j2=(j-1)*n+n;
                Sj=p(:,j1:j2)*p(:,j1:j2)'+eps0*eye(nd);
                scal=(det(Sj))^(1/nd); 
                S{j}=scal*inv(Sj);
            end   
        end
 
       % Compute auxiliares variables : D,P
       xv=repmat(x,1,c)-v(:,ic');
       ux=u(:)'; ux=ux(ones(nd,1),:); ux2=ux.*ux;
       bd=ux.*lam2-lam1+r*xv;
       add=r+r*ux2;
       bp=-lam2+r*bd.*ux./add;
       rux=r-r*r*ux2./add;
       for j=1:c
           j1=(j-1)*n+1; j2=(j-1)*n+n;
           ruxj=rux(:,j1:j2);
           bpj=bp(:,j1:j2);
           App=kron(speye(n),2*S{j})+spdiags(ruxj(:),0,nd*n,nd*n);
           ppx=App\bpj(:);
           p(:,j1:j2)=reshape(ppx,nd,n);
       end
       bd=r*ux.*p+bd;
       d=bd./add;
    
    end
    
    % Update multipliers
    lam1=lam1+r*(d-xv);
    lam2=lam2+r*(p-ux.*d);   

    % Compute errors
    ndu=sum(sum((u-uu).^2));   nu=sum(sum(u.^2));
    ndv=sum(sum((v-vv).^2));   nv=sum(sum(v.^2));
    nordd=sum(sum((d-dd).^2)); nord=sum(sum(d.^2)); 
    nordp=sum(sum((p-pp).^2)); norp=sum(sum(p.^2)); 
    nS=0; ndS=0;
    for j=1:c
        nS=nS+norm(S{j},'fro')^2;
        ndS=ndS+norm(S{j}-SS{j},'fro')^2;
    end 
    errd=[ndu;ndv;nordp;nordd;ndS];
    errn=[nu;nv;nord;norp;nS];
    
    err=sqrt(sum(errd)/sum(errn));

    if (iprint_inside == 1)    
        %Objectif function
        fobj = 0;
        for j=1:c
           Sj = S{j};
           for i=1:n
              pij = p(:,(j-1)*n+i); dij = d(:,(j-1)*n+i);
              fobj = fobj + pij'*Sj*pij;
           end
        end
       fprintf(">>iter=%i | err=%1.6e\n | J_FCM=%15.8e\n",iter,err,fobj);
       fprintf("ARI  = %10.4f \n",ARI(HP,Fuzzy2Hard(u)));

    end
 
end

%Exit verification : verification the ADMM problem converge to the FCM one.
max_residu = 10^-10;%maximal value of the residu
for j=1:c
   for i=1:n
       dij=d(:,(j-1)*n+i);pij=p(:,(j-1)*n+i);
       max_residu = max(max(max_residu,norm(pij-u(i,j)*dij)),norm(dij-x(:,i)+v(:,j)));
   end
end
   
if (iprint == 1) 

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
    fprintf(" J_FCM=%1.3e\n",fobj);
    fprintf("[iter=%i (max. %i)| err=%1.6e]\n",iter,itmax,err);    
    fprintf("Max residual value : %e \n",max_residu);
    if max_residu>10^-4
          fprintf("W - A - R - N - I - N - G : RESIDUAL VALUE TOO HIGH \n")
    end
    
    fprintf('*******************************************\n');
   
end
end

