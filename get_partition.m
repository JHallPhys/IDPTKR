function D1,get_partition(Hus_Entropy)

grid_size=size(Hus_Entropy);
Lq=1:grid_size(1); % Indices of the q(j) in [0,1)
Lp=Lq; % Same for p(j)
[Lqmesh,Lpmesh]=meshgrid(Lq,Lp); % Array of inidices of the elements grid

% Define the size of the box_size>=1

box_size=5;
dgrid=grid_size/box_size; % This is the number of points in grid that fit in the box in 1D 
bq=0:dgrid:grid_size 
bp=bq
[Bqmesh,Bpmesh]=meshgrid(bq,bp);

Partition=zeros(grid_size(1),grid_size(2));
for ittp=1:length(bp)-1
    for ittq=1:length(bq)-1
        Partition(:,:)=0;
        Left=Bqmesh(ittp,ittq)
        Right=Bqmesh(ittp,ittq+1)
        Up=Bpmesh(ittp,ittq);
        Down=Bpmesh(ittp+1,ittq);
        [i1,i2]=find(Lqmesh>Left & Lqmesh<Right & Lpmesh>Up & Lpmesh<Down);
        Partition(i1,i2)=1;
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
