function [] = DisplayClustering2D (X,V,HP,S,titre)
% 2D data (hard)partitioning display
%       AfficheClustering2D (X,V,HP,S,titre)
%
% INPUTS
%   X     : data (objects x attributs)
%   V     : centroïds (clusters x attributs)
%   HP    : Hard Partition (objects x 1)
%   S     : Mahalanobis matrix (cell clusters (attributs x attributs) ) 
%   Titre : Title de la figure (str)
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 03-25-2023
%  version: 1
%  --------------------------------------------------------------------------
  
  n = size(X,1); %Nb of objects
  c = size(V,1); %Nb of clusters
  p = size(V,2); %Nb of attributs
  
  %Projection : X,V -> 2D
  if p>2
    [axe1,axe2,pourcentage_info1,pourcentage_info2] = PCA_2D(X);
  else
      axe1 = [1; 0];axe2 = [0; 1];
      pourcentage_info1 = 100;pourcentage_info2 =100;
  end 
  XX = axe1'*X';YY = axe2'*X'; 
  Vx= axe1'*V';Vy= axe2'*V';
  
  %Unit circle (for S projections)
    n_u = 1000;
    t = linspace(0,2*pi,n_u);
    if mod(p,2)==0
        x = sin(t)/sqrt(p/2);y = cos(t)/sqrt(p/2);
        Cercle = [x' y'];
    else
        x = sin(t)/sqrt((p-1)/2); y = cos(t)/sqrt((p-1)/2);
        Cercle = [x' y' x'];
    end
    for k=2:p/2
        Cercle = [ Cercle , [x' y']];
    end
  
  %Assignement 
  XX_c = cell(c,1);YY_c = cell(c,1);
  for i=1:n
    XX_c{HP(i)}=[XX_c{HP(i)};XX(i)];
    YY_c{HP(i)}=[YY_c{HP(i)};YY(i)];
  end
   
  %Colours
  if c < 14 
     L_color = [[0,0,1];[1,0,0];[0,1,0];[0,1,1];[1,1,0];[1,0,1];...
      [0, 0.4470, 0.7410];[0.8500, 0.3250, 0.098];[0.929, 0.694, 0.125];...
      [0.494, 0.184, 0.556];[0.466, 0.674, 0.188];[0.301, 0.745, 0.933];...          	
      [0.635, 0.078, 0.184];[0, 0.5, 0];[0, 0.75, 0.75];[0.75, 0, 0.75];...
      [0.75, 0.75, 0];[0.25, 0.25, 0.25]];
  else
      L_color = [];for j=1:c;L_color = [L_color;[rand,rand,rand]]; end
  end
  
  %Starting - Display
  figure()
  
  labels={};
  hold on
  
  %For every cluster
  for j=1:c
      
    Xj=XX_c{j}; Yj=YY_c{j};  Vj = [Vx(j);Vy(j)];
 
    %Projection : S -> 2D
    [Q,D]=eig(S{j});        
    for i=1:n_u
      Ell_x(i) = axe1'*(Q*inv(sqrt(D))*Cercle(i,:)')+Vj(1);
      Ell_y(i) = axe2'*(Q*inv(sqrt(D))*Cercle(i,:)')+Vj(2);
    end
    
    %Test if the cluster has objects
    if max(size(Yj))~=0
        pX = plot(Xj,Yj,'o','MarkerSize',7); pX.Color = L_color(j,:);
        pV = plot(Vx(j),Vy(j),'o','MarkerSize',10);pV.MarkerFaceColor = L_color(j,:);
        pEll = plot(Ell_x,Ell_y,'--');pEll.Color = L_color(j,:);
        labels={labels{:},strcat('$\mbox{\boldmath{$\omega$}}_',num2str(j),'$'),...
        strcat('$\mbox{\boldmath{$v$}}_',num2str(j),'$'),strcat('$S_',num2str(j),'$')};      
    else
        disp(strcat('No assignement in the cluster',num2str(j)));
        pV = plot(Vx,Vy,'o');pV.Color = 'k';
        pV.MarkerFaceColor = L_color(j,:);pEll.Color = L_color(j,:);
        labels={labels{:},...
        strcat('$\mbox{\boldmath{$v$}}_',num2str(j),'$'),strcat('$S_',num2str(j),'$')};
    end
      
  end
  
  legend(labels,'Location','southwest','Interpreter','latex');
  xlabel(strcat('-- [',num2str(pourcentage_info1,'%.0f'),'%] -->'));
  ylabel(strcat('-- [',num2str(pourcentage_info2,'%.0f'),'%] -->'));
  title(titre);
  hold off
end