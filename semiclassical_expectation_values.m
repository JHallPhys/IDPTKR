clear all 
close all

n_efn=1;
%==========================================================================
%  Get the single state density
%==========================================================================

% User parameters
NC=1001;
k=2*pi;
gamma=0.001;
t_final=20;
hbar_half_sqrt=sqrt(1/(4*pi*NC));
sigma=NC*hbar_half_sqrt;
norm_index='FWD'
% norm_index='BWD';

% Initialisation
q = linspace(0,1,NC); 
p = linspace(-0.5,0.5,NC); 
dq=abs(q(2)-q(1));
[qmesh,pmesh]=meshgrid(q,p);


% Construct Norm map
[Norm_hm,Norm_hm_av] = nmap(t_final,qmesh,pmesh,NC,k,gamma,'BWD');

Norm_sort=Norm_hm_av(:);
Norm_sort=sort(Norm_sort,'descend');

CD=zeros(NC,NC);



for itt_efn=1:n_efn
itt_efn
tic 


if itt_efn==1

    Norm_single_state=Norm_hm_av;
    Norm_single_state(Norm_single_state<=Norm_sort(itt_efn*NC))=0;
    Norm_single_state(Norm_single_state>1)=1;
    CD_single=imgaussfilt(Norm_single_state,sigma);
    
else

    Norm_single_state=Norm_hm_av;
    Norm_single_state(Norm_single_state<=Norm_sort(itt_efn*NC))=0;
    Norm_single_state(Norm_single_state>Norm_sort((itt_efn-1)*NC))=0;
    Norm_single_state(Norm_single_state>1)=1;
    CD_single=imgaussfilt(Norm_single_state,sigma);

end

sum(sum(pmesh.*CD_single*dq^2))/gamma
    figure(1)
    clf
    imagesc(q,p,CD_single)
    % colorbar
    % caxis([0 1])
    colormap(viridis)
%     title(strcat('state',num2str(itt)))
    set(gca,'YDir','normal')


end

