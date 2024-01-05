clear all 
close all

n_efn=537
dim=2
ystr=strcat('D_',num2str(dim))
for itt_efn=1:n_efn
itt_efn
tic 
fname_efn=strcat('Husimi_Entropy_k10_g0p001_N1201_single_efn',num2str(itt_efn),'_special');

%==========================================================================
%  Load the quantum data
%==========================================================================
N=21;
parent_d = cd;    
cd './Husimi_dat' % Directory where matrix is stored
Hus_Entropy = matfile(fname_efn);
Hus_Entropy=Hus_Entropy.Hus_Entropy; % I think this step may be redundent
cd(parent_d)
% figure
% imagesc(Hus_Entropy)
% colorbar
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% % caxis([0 1])
% c = colorbar('eastoutside');
% Hus_Entropy=Hus_Entropy;
% return
% sum(sum(Hus_Entropy))
% return
%==========================================================================
%  Get the partition of the phase space 
%==========================================================================
% 
% You can hardcode Lq,Lp directly into meshgrid to cut all this to 2 lines
grid_size=size(Hus_Entropy);
grid_size=grid_size(1);
% Hus_Entropy=Hus_Entropy./grid_size(1);
Lq=1:grid_size; % Indices of the q(j) in [0,1)
Lp=Lq; % Same for p(j)
[Lqmesh,Lpmesh]=meshgrid(Lq,Lp); % Array of inidices of the elements grid

% Define the size of the box_size>=0
% Nmin=1; %
% Nmax=round(sqrt(grid_size)); % Maximum box size in 1d <sqrt(1/h)=sqrt(N)

Nmin=1; %
Nmax=sqrt(grid_size);% Maximum box size in 1d <sqrt(1/h)=sqrt(N)
N_i=linspace(Nmin,Nmax,N);
% N_i=fliplr(linspace(Nmin,Nmax,N));
dN=abs(N_i(2)-N_i(1)); % This is the spacing in the boxes
NLen=length(N_i);
% N_i(1)=N_i(1)+1
% return

SE=zeros(NLen,1); % Entropy Array
Nmax=NLen-1;
% Nmax=10
% for itt_box=1:NLen


%==========================================================================
%==========================================================================

for itt_box=1:Nmax
    itt_box;

    box_size=N_i(itt_box); % Size of the box (1/N_i)^2
    dgrid=grid_size/box_size; % This is the number of points in grid that fit in the box in 1D 
    bq=0:dgrid:grid_size; % Boxes placed on the interval 
    bp=bq;
    [Bqmesh,Bpmesh]=meshgrid(bq,bp);
    Partition=zeros(grid_size,grid_size);
   

%     return
    % Check if final element of bq is equal to grid_size
    eps_grid=1e-6; % Parameter to check above condition
 

    %==========================================================================
    % Grid is fine
    %==========================================================================


    if abs(grid_size-bq(end))<eps_grid % Pass
    
        SE=box_measure_higher(SE,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy,dgrid,dim);

    end
    
    %==========================================================================
    % Need Bigger Grid
    %==========================================================================
    
    
    if abs(grid_size-bq(end))>eps_grid % Pass
        
        bq(length(bq)+1)=bq(length(bq))+dgrid; % Add a new element outside of the grid
        bp=bq; % Same for this one 
        [Bqmesh,Bpmesh]=meshgrid(bq,bp); % Make new meshgrid
%         return
        SE=box_measure_higher(SE,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy,dgrid,dim);
%     return
    end
    
    SE(itt_box,1);
% return


end


% figure
% plot((log(1./N_i(1:Nmax))),SE(1:Nmax),'k.-','markersize',10)
% xlabel('\epsilon')
% ylabel('S_1(\epsilon)')


DC = polyfit(log(1./N_i(1:Nmax)), SE(1:Nmax), 1);
D(itt_efn)=DC(1);

toc

figure(1)
clf
plot(1:1:itt_efn,D,'k.-','markersize',5)
xlabel('\gamma')
ylabel(ystr)
axis([1 n_efn 1.5 2])

figure(2)
clf
plot(1:1:itt_efn,D,'k.','markersize',10)
xlabel('\gamma')
ylabel(ystr)
axis([1 n_efn 1.5 2])

end


figure
clf
plot(1:1:itt_efn,D,'k.-','markersize',10)
xlabel('\gamma')
ylabel(ystr)
axis([1 itt_efn 1.5 2])


figure
clf
plot(1:1:itt_efn,D,'k.','markersize',10)
xlabel('n_{\gamma}')
ylabel(ystr)
axis([1 itt_efn 1.5 2])




