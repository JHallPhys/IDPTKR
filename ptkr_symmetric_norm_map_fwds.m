clear all
close all

% User parameters
N=1001;
k=5;
gamma=0.001;
t_final=18;

norm_index='FWD'
% norm_index='BWD';

% Initialisation
q = linspace(0,1,N); 
p = linspace(-0.5,0.5,N); 
[qmesh,pmesh]=meshgrid(q,p); 

% Construct Norm map
[Norm_hm,Norm_hm_av] = nmap(t_final,qmesh,pmesh,N,k,gamma,norm_index);
% Norm_hm=Norm_hm_av;
% 
inferno=inferno();
figure
imagesc(q,p,Norm_hm)
colorbar
colormap(inferno)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
% caxis([0 1])
c = colorbar('eastoutside');
