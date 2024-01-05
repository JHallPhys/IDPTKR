clear all 
close all

%==========================================================================
%  Load the quantum data
%==========================================================================
N=101
parent_d = cd;    
cd './Husimi_dat' % Directory where matrix is stored
Hus_Entropy = matfile('Husimi_Entropy_k10_g0p001_N2001_single_efn5');
Hus_Entropy=Hus_Entropy.Hus_Entropy; % I think this step may be redundent
cd(parent_d)


%==========================================================================
%  Get the partition of the phase space 
%==========================================================================
% 
% You can hardcode Lq,Lp directly into meshgrid to cut all this to 2 lines
grid_size=size(Hus_Entropy);
grid_size=grid_size(1);
% Hus_Entropy=Hus_Entropy./grid_size(1);
Lq=1:grid_size(1); % Indices of the q(j) in [0,1)
Lp=Lq; % Same for p(j)
[Lqmesh,Lpmesh]=meshgrid(Lq,Lp); % Array of inidices of the elements grid

% Define the size of the box_size>=0
Nmin=4; %
Nmax=round(sqrt(grid_size(1))); % Maximum box size in 1d <sqrt(1/h)=sqrt(N)
N_i=linspace(Nmin,Nmax,N);
dN=abs(N_i(2)-N_i(1)) % This is the spacing in the boxes
NLen=length(N_i);
% N_i(1)=N_i(1)+1
% return

SE=zeros(NLen,1); % Entropy Array
Nmax=NLen-1;
Nmax=100
% for itt_box=1:NLen
for itt_box=2:Nmax
    itt_box

    box_size=N_i(itt_box); % Size of the box (1/N_i)^2
    dgrid=grid_size/box_size; % This is the number of points in grid that fit in the box in 1D 
    bq=0:dgrid:grid_size; % Boxes placed on the interval 
    bp=bq;
    [Bqmesh,Bpmesh]=meshgrid(bq,bp);
    Partition=zeros(grid_size,grid_size);
    tic

    
    % Check if final element of bq is equal to grid_size
    eps_grid=1e-6; % Parameter to check above condition
 

    %==========================================================================
    % Pass
    %==========================================================================


    if abs(grid_size-bq(end))<eps_grid % Pass
    
        SE=box_grid(SE,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy);

    end
    
    %==========================================================================
    % Fail
    %==========================================================================
    
    
    if abs(grid_size-bq(end))>eps_grid % Pass
        
        bq(length(bq)+1)=bq(length(bq))+dgrid; % Add a new element outside of the grid
        bp=bq; % Same for this one 
        [Bqmesh,Bpmesh]=meshgrid(bq,bp); % Make new meshgrid
%         return
        SE=box_grid(SE,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy);
%     return
    end



end

                                                    

