clear all 
close all


parent_d = cd;    
cd './Husimi_dat' % Directory where matrix is stored
% Hus_Entropy = matfile('Husimi_Entropy_k10_g0p001_N1001_subset60');
% Hus_Entropy = matfile('Husimi_Entropy_k5_g0p001_N1001_subset40');
% Hus_Entropy = matfile('Husimi_Entropy_k10_g0p001_N1001_subset100');
% Hus_Entropy = matfile('Husimi_Entropy_k10_g0p001_N2001_subset20');
Hus_Entropy = matfile('Husimi_Entropy_k10_g0p001_N2001_single_efn5');
% Hus_Entropy=matfile('Husimi_Entropy_k10_g0p001_N2001_subset60');
% Hus_Entropy=matfile('Husimi_Entropy_k10_g0p001_N1001_subset5');
% Hus_Entropy=matfile('Husimi_Entropy_k10_g0p001_N1001_subset3');
Hus_Entropy=Hus_Entropy.Hus_Entropy; % I think this step may be redundent
cd(parent_d)

% Normp=get_circle(0.3,qmesh,pmesh,N);



% figure
% imagesc(Hus_Entropy)
% colorbar
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% % caxis([0 1])
% c = colorbar('eastoutside');

% return

%==========================================================================
%  Get the partition of the phase space 
%==========================================================================
% 
% You can hardcode Lq,Lp directly into meshgrid to cut all this to 2 lines
grid_size=size(Hus_Entropy);
% Hus_Entropy=Hus_Entropy./grid_size(1);
Lq=1:grid_size(1); % Indices of the q(j) in [0,1)
Lp=Lq; % Same for p(j)
[Lqmesh,Lpmesh]=meshgrid(Lq,Lp); % Array of inidices of the elements grid

% Define the size of the box_size>=1
Nmin=0; %
dN=0.5;
Nmax=round(sqrt(grid_size(1))) 
N_i=Nmin:dN:Nmax;
NLen=length(N_i);
% N_i(1)=N_i(1)+1
return

SE=zeros(NLen,1) % Entropy Array

% for itt_box=1:NLen
for itt_box=1:NLen

box_size=N_i(itt_box)
dgrid=grid_size/box_size; % This is the number of points in grid that fit in the box in 1D 
bq=0:dgrid:grid_size ;
bp=bq;
[Bqmesh,Bpmesh]=meshgrid(bq,bp);
Partition=zeros(grid_size(1),grid_size(2));
tic
for ittp=1:length(bp)-1
    for ittq=1:length(bq)-1
%         Partition(:,:)=0;
        Left=Bqmesh(ittp,ittq);
        Right=Bqmesh(ittp,ittq+1);
        Up=Bpmesh(ittp,ittq);
        Down=Bpmesh(ittp+1,ittq);
        [i1,i2]=find(Lqmesh>Left & Lqmesh<Right & Lpmesh>Up & Lpmesh<Down);
        Partition(min(i1):max(i1),min(i2):max(i2))=1;
        SE(itt_box,1)= SE(itt_box,1)+sum(sum(Hus_Entropy(min(i1):max(i1),min(i2):max(i2)))); % This should be multiplied by a measure but is it dqdp or deps^2?
        figure(2)
        clf
        imagesc(bq,bp,Partition)
        colorbar
        colormap(viridis)
        set(gca,'YDir','normal')
        xlabel('q')
        ylabel('p')
        % caxis([0 1])
        c = colorbar('eastoutside');
%         return
    end
end
return
toc
end

for itt=1:NLen-1
    D(itt)=(SE(itt)-SE(itt+1))/(log(N_i(itt))-log(N_i(itt+1)));
end

figure
plot(N_i(1:NLen-1),D,'k.')
return

figure
plot(log10(1./N_i),SE,'k.-','markersize',10)
xlabel('\epsilon')
ylabel('N(\epsilon)')

figure
plot(log(1./N_i),SE,'k.-','markersize',10)
xlabel('\epsilon')
ylabel('S_1(\epsilon)')


return

