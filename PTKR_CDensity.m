clear all 
close all
%==========================================================================
%  Get the single state density
%==========================================================================

% User parameters
N=1001;
k=10;
gamma=0.001;
hbar_half_sqrt=sqrt(1/(4*pi*N));
sigma=N*hbar_half_sqrt;
norm_index='FWD'
% norm_index='BWD';

% Initialisation
q = linspace(0,1,N); 
p = linspace(-0.5,0.5,N); 
dq=abs(q(2)-q(1));
[qmesh,pmesh]=meshgrid(q,p);


if N==501
    nfq=154/N
    t_final=40
end
if N==1001
    nfq=431/N
    t_final=40
end
if N==1501
    nfq=681/N
    t_final=10
%     return
end
if N==2001
    nfq=944/N
end

Hus=zeros(N,N);
tic
for itt_efn = 1:round(nfq*N)
    itt_efn
fname_efn=strcat('Husimi_Entropy_k5_g0p001_N',num2str(N),'_single_efn',num2str(itt_efn),'_special');
parent_d = cd;    
cd './Husimi_dat' % Directory where matrix is stored
Hus_Entropy = matfile(fname_efn);
Hus_Entropy=Hus_Entropy.Hus_Entropy; % I think this step may be redundent
cd(parent_d)
Hus=Hus+Hus_Entropy;
end
toc

    figure(1)
    clf
    imagesc(q,p,Hus)
    % colorbar
    caxis([0 1])
    colormap(viridis)
%     title(strcat('state',num2str(itt)))
    set(gca,'YDir','normal')


return
% Construct Norm map
[Norm_hm,Norm_hm_av] = nmap(t_final,qmesh,pmesh,N,k,gamma,'BWD');
Norm_sort=Norm_hm_av(:);
Norm_sort=sort(Norm_sort,'descend');

CD=zeros(N,N);


Norm_single_state=Norm_hm_av;
Norm_single_state(Norm_single_state<=Norm_sort(round(nfq*N)*N))=0;
Norm_single_state(Norm_single_state>1)=1;
CD_single=imgaussfilt(Norm_single_state,sigma);




sum(sum(pmesh.*CD_single*dq^2))/gamma;
    figure(1)
    clf
    imagesc(q,p,CD_single)
    % colorbar
    caxis([0 1])
    colormap(viridis)
%     title(strcat('state',num2str(itt)))
    set(gca,'YDir','normal')



FDG_nrm=FDG./(N^2*nfc);
M=FDG_nrm+Hus_nrm;
M=0.5*M;
Jdiv=0.5*(kldiv(FDG_nrm,M))+0.5*(kldiv(Hus_nrm,M));
% Entropy=kldiv(FDG_nrm,Hus_nrm);

% tfit(count,1)=Entropy;    
tfit(count,1)=Jdiv; 
tfit(count,2)=itt_time;



