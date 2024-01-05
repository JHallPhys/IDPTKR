function [norm_out,qmesh_out,pmesh_out] = normfwd(k,gamma,qmesh_in,pmesh_in)


% Partially frozen
%==========================================================================

% norm_out = exp(2*gamma*pmesh_in);
pstar_in=pmesh_in+k*sin(2*pi*qmesh_in)/(4*pi);
pnorm=mod(pstar_in+1/2,1)-1/2;
norm_out = exp(2*gamma*pnorm); % 
qmesh_in=mod(qmesh_in+pstar_in+gamma/(2*pi),1);
pmesh_in=mod(pstar_in+k*sin(2*pi*qmesh_in)/(4*pi)+gamma/(2*pi)+1/2,1)-1/2;    


    

qmesh_out = qmesh_in;
pmesh_out= pmesh_in;

end

