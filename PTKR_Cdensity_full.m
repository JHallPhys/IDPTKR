clear all 
close all
%==========================================================================
%  Get the single state density
%==========================================================================

% User parameters
N=3201;
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
    t_final=200
end
if N==1001
    nfq=431/N
    t_final=60
end
if N==1501
    nfq=681/N
    t_final=40
%     return
end
if N==2001
    nfq=944/N
    t_final=40
end
if N==2401
    nfq=1139/N
    t_final=40
end
if N==3201
    nfq=1539/N
    t_final=40
end


%==========================================================================
% Cdensity Main Algor
%==========================================================================

% Get Husimi and Normalise it

Hus=gain_set_check(k,N,gamma,nfq); % Load the gain set from the single state data
HN=sum(sum(Hus.*dq^2));
Hus_nrm=Hus./(HN);


% Calculate the JS Divergence for each time step

tfit=zeros(t_final,2); % array whose columns will be JS divergence and time

for itt_time=1:t_final

itt_time

[Norm_hm,Norm_hm_av] = nmap(itt_time,qmesh,pmesh,N,k,gamma,'BWD'); % Get norm landscape
CD=get_cd_ptkr(Norm_hm_av,sigma,nfq,N); % Get classical density
CN=sum(sum(CD.*dq^2)); % Normalisation constant for classical density
CD_nrm=CD./CN; % Normalised classical density
M=CD_nrm+Hus_nrm; % Distribution for Jdov
M=0.5*M;
Jdiv=0.5*(kldiv(CD_nrm,M))+0.5*(kldiv(Hus_nrm,M)); % Calculate JDiv  
tfit(itt_time,1)=itt_time; % Time
tfit(itt_time,2)=Jdiv*dq^2; % Jensen-Shannon divergence
% return

end

figure
plot(tfit(:,1),tfit(:,2),'k.-','markersize',2)
xlabel('t_f')
ylabel('I_{JS}')
axis([0 t_final 0 0.25])
set(gca, 'YScale', 'log')


[tbest,ind_best]=min(tfit(:,2)); % Find the best time

% Plot classical density at best time

[Norm_hm,Norm_hm_av] = nmap(tfit(ind_best,1),qmesh,pmesh,N,k,gamma,'BWD'); % Get norm landscape
CD_Best=get_cd_ptkr(Norm_hm_av,sigma,nfq,N); % Get classical density

figure
clf
imagesc(q,p,Hus)
colorbar
caxis([0 1])
colormap(viridis)
set(gca,'YDir','normal')
title('Quantum')


figure
clf
imagesc(q,p,CD_Best)
colorbar
caxis([0 1])
colormap(viridis)
set(gca,'YDir','normal')
title('Classical')

'best time'
tfit(ind_best,1)