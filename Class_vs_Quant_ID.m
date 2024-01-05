clear all 
close all

n_efn=400;
%==========================================================================
%  Get the single state density
%==========================================================================

% User parameters
NC=1001;
k=5;
gamma=0.001;
t_final=60;
hbar_half_sqrt=sqrt(1/(4*pi*NC));
sigma=NC*hbar_half_sqrt;
norm_index='FWD'
% norm_index='BWD';

% Initialisation
q = linspace(0,1,NC); 
p = linspace(-0.5,0.5,NC); 
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

fname_efn=strcat('Husimi_Entropy_k5_g0p001_N1001_single_efn',num2str(itt_efn),'_special');
parent_d = cd;    
cd './Husimi_dat' % Directory where matrix is stored
Hus_Entropy = matfile(fname_efn);
Hus_Entropy=Hus_Entropy.Hus_Entropy; % I think this step may be redundent
cd(parent_d)

%==========================================================================
%  Load the quantum data
%==========================================================================
N=21;
Hus_Entropy_class=CD_single;
Hus_Entropy_quant=Hus_Entropy;
%==========================================================================
%  Get the partition of the phase space 
%==========================================================================
% 
% You can hardcode Lq,Lp directly into meshgrid to cut all this to 2 lines
grid_size=size(Hus_Entropy_class);
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

SE_class=zeros(NLen,1); % Entropy Array
SE_quant=zeros(NLen,1); % Entropy Array
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
    
        SE_class=box_measure(SE_class,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy_class,dgrid);
        SE_quant=box_measure(SE_quant,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy_quant,dgrid);

    end
    
    %==========================================================================
    % Need Bigger Grid
    %==========================================================================
    
    
    if abs(grid_size-bq(end))>eps_grid % Pass
        
        bq(length(bq)+1)=bq(length(bq))+dgrid; % Add a new element outside of the grid
        bp=bq; % Same for this one 
        [Bqmesh,Bpmesh]=meshgrid(bq,bp); % Make new meshgrid
%         return
        SE_class=box_measure(SE_class,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy_class,dgrid);
        SE_quant=box_measure(SE_quant,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy_quant,dgrid);
%     return
    end
    
    SE_class(itt_box,1);
% return


end


% figure
% plot((log(1./N_i(1:Nmax))),SE(1:Nmax),'k.-','markersize',10)
% xlabel('\epsilon')
% ylabel('S_1(\epsilon)')



DC_class = polyfit(log(1./N_i(1:Nmax)), SE_class(1:Nmax), 1);
IDClass(itt_efn)=DC_class(1);


DC_quant = polyfit(log(1./N_i(1:Nmax)), SE_quant(1:Nmax), 1);
IDQuant(itt_efn)=DC_quant(1);

toc

figure(1)
clf
hold on
plot(1:1:itt_efn,IDClass,'r.','markersize',10)
plot(1:1:itt_efn,IDQuant,'b.','markersize',10)
xlabel('\gamma')
ylabel('D_1(\epsilon)')
axis([1 n_efn 1.5 2])

figure(2)
clf
hold on
plot(1:1:itt_efn,IDClass,'r.-','markersize',5)
plot(1:1:itt_efn,IDQuant,'b.-','markersize',5)
xlabel('\gamma')
ylabel('D_1(\epsilon)')
axis([1 n_efn 1.5 2])

end







