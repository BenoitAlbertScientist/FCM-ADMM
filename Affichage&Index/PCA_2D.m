function [axe1,axe2,pourcentage_info1,pourcentage_info2] = PCA_2D(X)
% Principal component analysis
%       [axe1,axe2,pourcentage_info1,pourcentage_info2] = PCA_2D(X)
%
% INPUTS
%   X : dataset (objects x  attributs)
%
% OUTPUTS
%   axe1
%   axe2 
%   pourcentage_info1
%   pourcentage_info1
%
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 10-21-2021
%  version: 2
%  --------------------------------------------------------------------------

    % Number of objects in the dataset
    n = size(X, 1); 
    
    % Centring the data by subtracting the mean of each dimension
    X_ = X - mean(X) .* ones(n, 1);
    
    % Calculation of eigenvectors (axis1 and axis2) and eigenvalues
    [Q, D] = eig(X_' * X_);
    
    % Extraction of eigenvalues in a vector
    D = diag(D)';
    
    % Selection of the two largest eigenvalues and their indices
    [vals, composante_principales] = maxk(D, 2);
    
    % Extraction of the first two principal components (axis1 and axis2)
    axe1 = Q(:, composante_principales(1));
    axe2 = Q(:, composante_principales(2));
    
    % Calculation of the percentage of information carried by each principal component
    somme_totale_vp = sum(D);
    pourcentage_info1 = (vals(1) / somme_totale_vp) * 100;
    pourcentage_info2 = (vals(2) / somme_totale_vp) * 100;
    
end

