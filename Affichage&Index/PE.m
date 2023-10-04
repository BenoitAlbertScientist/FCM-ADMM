function [pe] = PE(U)
% Partition Entropy
% https://en.wikipedia.org/wiki/Rand_index
%        score_ari = PE(U)
%
% INPUTS
%   U  : Fuzzy partition - probability  (clusters x objects)
% OUTPUTS
%   pe : partition entropy
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 10-01-2023
%  version: 0
%  --------------------------------------------------------------------------

    [n,c]=size(U);
    pe = 0; 
    U=max(U,0.000000001*ones(n,c)); %Avoid U(i,j) = 0 => log(U(i,j)) = NAN

    for j=1:c 
        pe = pe - U(:,j)'*log2(U(:,j));
    end
    pe = pe/n;

end

