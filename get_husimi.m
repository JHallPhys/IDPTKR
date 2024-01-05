function Hus=get_husimi(N,n_efn,q,p,z,psi_gamma)

Hus = zeros(N,N,n_efn);
cs=zeros(N,N);
norm_cs= (2/N)^0.25; % Normalisation constant for the coherent state
tic
% return
for itt = 1:N-1
    itt
   
    cs=Cs_create_component(itt,norm_cs,N,q,z,cs); 
    

    
    phi2(1,1,:)=psi_gamma(itt,:);
    Hus=Hus+conj(cs).*phi2;
    

    cs(:,:)=0;
end
toc






end