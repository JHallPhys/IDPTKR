function [y1,y2,y3] = normbwd(k,gamma,qmesh,pmesh)


% Partially frozen
%==========================================================================
% y1 = exp(-2*gamma*pmesh);
% y1 = exp(-2*gamma*pmesh);


pstar=pmesh-k*sin(2*pi*qmesh)/(4*pi);
pnorm=mod(pstar+1/2,1)-1/2;
norm_out = exp(2*gamma*pnorm); %
qmesh=mod(qmesh-pstar+gamma/(2*pi),1);
pmesh=mod(pstar-k*sin(2*pi*qmesh)/(4*pi)+gamma/(2*pi)+1/2,1)-1/2;    
   


y1=norm_out;
y2 = qmesh;
y3= pmesh;

end

