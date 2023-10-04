function [score_ari,score_ri] = ARI (HP1, HP2)
% Ajusted Rand Index between two (Hard) partitions
% https://en.wikipedia.org/wiki/Rand_index
%        score_ari = ARI (Hp1, Hp2)
%
% INPUTS
%   HP1, HP2 : Hard partitions - (objects x 1)
% OUTPUTS
%   score_ari
%   score_ri
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 09-25-2023
%  version: 2
%  --------------------------------------------------------------------------

n=length(HP1);
n2=length(HP2);
if n~=n2
    disp('No same number of objects ')
end

a=0;b=0;c=0;d=0;

for i=1:n
  for j=i+1:n
    if HP2(i)==HP2(j)
      if HP1(i)==HP1(j)
        a=a+1;
      else
        b=b+1;
      end
    else
      if HP1(i)==HP1(j)
        c=c+1;
      else
        d=d+1;
      end
    end
  end
end

Index = a;
ExpectedIndex = ((a+b)*(a+c))/(a+b+c+d);
MaxIndex = ((a+b)+(a+c))/2;

score_ari = (Index - ExpectedIndex)/(MaxIndex - ExpectedIndex);
score_ri = (a+d)/(a+b+c+d);

end

