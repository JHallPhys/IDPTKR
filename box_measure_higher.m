function SE=box_measure_higher(SE,itt_box,Partition,bp,bq,Bpmesh,Bqmesh,Lpmesh,Lqmesh,Hus_Entropy,dgrid,dim)
  L=length(Partition);
  Husimi_partition=zeros(L,L);
  dHus=1/sum(sum(Hus_Entropy));
  
    for ittp=1:length(bp)-1
        for ittq=1:length(bq)-1

            % Define the corners of the box
            Left=Bqmesh(ittp,ittq);
            Right=Bqmesh(ittp,ittq+1);
            Up=Bpmesh(ittp,ittq);
            Down=Bpmesh(ittp+1,ittq);

            % Finsd the elements inside the box
            [i1,i2]=find(Lqmesh>Left & Lqmesh<Right & Lpmesh>Up & Lpmesh<Down);
%             Hus_Entropy(min(i1):max(i1),min(i2):max(i2))=1;
            Pbox=sum(sum(Hus_Entropy(min(i1):max(i1),min(i2):max(i2))))*dHus;
%              Pbox=sum(sum(Hus_Entropy(min(i1):max(i1),min(i2):max(i2)).*log(Hus_Entropy(min(i1):max(i1),min(i2):max(i2)))))*dz; % The Husimi measure inside the partition
%              Pbox=Pbox; % THe probability inside the partition
             if Pbox<1e-12
                  SE(itt_box,1)= SE(itt_box,1);
             else
                SE(itt_box,1)= SE(itt_box,1)+Pbox^dim; % This should be multiplied by a measure but is it dqdp or deps^2?
             end
%            SE(itt_box,1)= SE(itt_box,1)-Pbox; % This should be multiplied by a measure but is it dqdp or deps^2?
%              SE(itt_box,1)= SE(itt_box,1)+Pbox*log(Pbox)*(dgrid); % This should be multiplied by a measure but is it dqdp or deps^2?

             Husimi_partition(min(i1):max(i1),min(i2):max(i2))= Husimi_partition(min(i1):max(i1),min(i2):max(i2))+Hus_Entropy(min(i1):max(i1),min(i2):max(i2));
            
%             return
%              SE(itt_box,1)= SE(itt_box,1)+sum(sum(Hus_Entropy(min(i1):max(i1),min(i2):max(i2))*dgrid)); % This should be multiplied by a measure but is it dqdp or deps^2?
%             figure(2)
%             clf
%             imagesc(bq,bp,Husimi_partition)
%             colorbar
%             colormap(viridis)
%             set(gca,'YDir','normal')
%             xlabel('q')
%             ylabel('p')
% %             caxis([0 1])
%             c = colorbar('eastoutside');
            
        
%             return
        end
    end
    
%     rerturn
SE(itt_box,1)=-log(SE(itt_box,1))./(1-dim);
end