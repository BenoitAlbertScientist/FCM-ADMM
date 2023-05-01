function [axe1,axe2] = PCA_2D(X)
%  Principal component analysis n dimensions -> 2D
%  https://en.wikipedia.org/wiki/Principal_component_analysis
%   [axe1,axe2] = PCA_2D(X)
%
% INPUTS
%   X : data (clusters x objects)
%
% OUTPUTS
%   axe1,axe2 : axes in 2D
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 03-25-2023
%  version: 1
%  --------------------------------------------------------------------------


   n = size(X,1);
   
   X_ = X - mean(X).*ones(n,1);
   
   [Q,D] = eig(X_'*X_);
   
   D = diag(D)';
   [vals,composante_principales] = maxk(D,2);
   
   axe1 = Q(:,composante_principales(1));
   axe2 = Q(:,composante_principales(2));
   
  
end

