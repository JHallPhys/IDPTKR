function [div_out] = kldiv(A_in,B_in)

% Calculates the Kulllback-Lieber Divergence for finite probability distributions
% A_in and B_in via the formula:
%
% Div=sum A_in log2(A_in/B_in)



N=length(A_in);
div_out=0;

for j1 = 1:N
for k2 = 1:N
    
    if A_in(j1,k2)<=1e-12
        div_out=div_out;
    else
      div_out=div_out+A_in(j1,k2)*log2(A_in(j1,k2)/B_in(j1,k2)); 
    end
end


end
