function [fs]=FS(X, U)
% Fuzzy Silhouette index
%   [fs]=FS(X, U)
%
% INPUTS
%   X  : dataset (objects x  attributs)
%   U  : Fuzzy partition (clusters x objects)
%
% OUTPUTS
%   fs : fussy silhouette index
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 03-25-2023
%  version: 1
%  --------------------------------------------------------------------------

alpha = 1; % Coefficiant of fuzziness
[n,c]=size(U);
HP = Fuzzy2Hard(U);

sil = zeros(n,1);u_u = zeros(n,1);
for i = 1:n
    
    %Select two highest membership degree and n° clusters associated
    [u_i,j]=maxk(U(i,:), 2);
    
    a = 0;b = 10^16;

    list_obj_wj1 = find(HP==j(1));%Element de la classe k
    if length(list_obj_wj1)>1 % If the cluster isn't empty
        som_d = 0;
        for ind_l=1:length(list_obj_wj1)
            l = list_obj_wj1(ind_l);
            som_d = som_d + sqrt((X(i,:)-X(l,:))*(X(i,:)-X(l,:))');
        end
     a = som_d/(length(list_obj_wj1)-1);
    end

    list_obj_wj2 = find(HP==j(2));%Element de la classe k
    if length(list_obj_wj2)>0 % If the cluster isn't empty
        som_d = 0;
        for ind_l=1:length(list_obj_wj2)
            l = list_obj_wj2(ind_l);
            som_d = som_d + sqrt((X(i,:)-X(l,:))*(X(i,:)-X(l,:))');
        end
        b = som_d/length(list_obj_wj2);    
    end
    sil(i) = (b-a)/max(b,a);
    u_u(i) = (u_i(1)-u_i(2))^alpha;
end

fs = sum(u_u.*sil)/sum(u_u);

end

function [HP] = Fuzzy2Hard(U)
% Transform a fuzzy partition (U) in to an hard partition (HP)

    HP = zeros(size(U,1),1);
    for k = 1:size(U,1)
        [vv,HP(k)] = max(U(k,:));
    end
    
end