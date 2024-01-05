function [Norm_hm_out,Norm_hm_av_out] = nmap(t_final,qmesh_in,pmesh_in,N_in,k_in,gamma_in,norm_index_in)

%==========================================================================
% Allocation
%==========================================================================

qmesh=qmesh_in;
pmesh=pmesh_in;
Norm_hm=zeros(N_in,N_in); 
Norm_hm_av = zeros(N_in,N_in);
Norm_hm_0=zeros(N_in,N_in);
Norm_hm_0(:,:)=1;
k=k_in;
gamma=gamma_in;
%==========================================================================
% Construction
%==========================================================================

if isequal(norm_index_in,'FWD') % Forward norm map
    
    for t = 1:t_final

    [Norm_hm,qmesh,pmesh] = normfwd(k,gamma,qmesh,pmesh);
    Norm_hm_0=Norm_hm.*Norm_hm_0;
    Norm_hm_av = Norm_hm_av + Norm_hm_0;

    end
    
    
elseif isequal(norm_index_in,'BWD') % Backward norm map

    for t = 1:t_final

    [Norm_hm,qmesh,pmesh] = normbwd(k,gamma,qmesh,pmesh);
    Norm_hm_0=Norm_hm.*Norm_hm_0;
    Norm_hm_av = Norm_hm_av + Norm_hm_0;
    
    end
    
end



Norm_hm_out=Norm_hm_0;
Norm_hm_av_out=Norm_hm_av./t_final;

end