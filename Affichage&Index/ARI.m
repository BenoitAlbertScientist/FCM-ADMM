function score_ari = ARI (HP1, HP2)
% Ajusted Rand Index between two (Hard) partitions
% https://en.wikipedia.org/wiki/Rand_index
%        score_ari = ARI (Hp1, Hp2)
%
% INPUTS
%   HP1, HP2 : Hard partition - array (n,1)
% OUTPUTS
%   score_ari
%
%  --------------------------------------------------------------------------
%  Author : Benoit Albert
%  mail   : benoit.albert@uca.fr
%  date   : 03-25-2023
%  version: 1
%  --------------------------------------------------------------------------

a=0;b=0;c=0;d=0;

n=length(HP2);
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
end
