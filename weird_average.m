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
    
    t_final=60
end
if N==1001
    nfq=431/N
    nfq=1/N;
    

    t_final=50

end
if N==1501
    nfq=681/N
    t_final=20
%     return
end
if N==2001
    nfq=944/N
    t_final=16
end

CD_av=zeros(N,N);

for itt=1:100
itt
% Construct Norm map
[Norm_hm,Norm_hm_av] = nmap(t_final,qmesh+rand(N,N)./N,pmesh+rand(N,N)./N,N,k,gamma,'BWD');
Norm_sort=Norm_hm_av(:);
Norm_sort=sort(Norm_sort,'descend');

CD=zeros(N,N);


Norm_single_state=Norm_hm_av;
Norm_single_state(Norm_single_state<=Norm_sort(round(nfq*N)*N))=0;
Norm_single_state(Norm_single_state>1)=1;
CD_single=imgaussfilt(Norm_single_state,sigma);

CD_av=CD_av+CD_single;

sum(sum(pmesh.*CD_single*dq^2))/gamma;
figure(1)
clf
imagesc(q,p,CD_av./itt)
colorbar
% caxis([0 1])
colormap(viridis)
%     title(strcat('state',num2str(itt)))
set(gca,'YDir','normal')

end





