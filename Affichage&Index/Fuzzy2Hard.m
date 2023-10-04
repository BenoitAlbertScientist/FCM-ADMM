function [HP] = Fuzzy2Hard(U)
% Transform a fuzzy partition (U) in to an hard partition (HP)
% Assignment to the cluster with the highest degree of membership
%   [HP] = Fuzzy2Hard(U) 
%
% INPUTS
%   U  : Fuzzy partition - probability  (clusters x objects)
%
% OUTPUTS
%   HP : Hard partition (objects x 1) 
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 03-25-2023
%  version: 1
%  --------------------------------------------------------------------------


    HP = zeros(size(U,1),1);
    for k = 1:size(U,1)
        [vv,HP(k)] = max(U(k,:));
    end
    
end