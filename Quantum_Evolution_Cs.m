clear all
close all
parallel.gpu.enableCUDAForwardCompatibility(true) % RTX 3090 ;)
%==========================================================================
%   Quantum System and Dynamical Parameters
%==========================================================================

% System Parameters

N_1 =250; 
N = 2*N_1+1; % Hilbert space dimensio
K_class =0.5; % Classical Kicking
gamma = complex(0,0.1); % PT-strength 
q_0 = 0.5; % 
p_0 = 0.2;
t_final=10;

T=2*pi/N; % Effective hbar
kick = K_class/T; % Quantum Kicking

% Dynamics Parameters 

Srep='n'%
set_efn='G'; % Invarient Subspace: Gain ('G'), Stable ('S'), Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Stable ('0'), Loss ('-')
% return
%==========================================================================
%   Matrix Construction and Schur decompesition
%==========================================================================

U=zeros(N,N); % Initialise Flouqet matrix
[U,time1]=UMatrix(U,N,N_1,K_class,T,gamma); % Construct Flouqet matrix
[psi,En] = schur(U); % psi are the Schur eigenfns and En matrix of eigs
[psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values

%==========================================================================
%  Set up the grid
%==========================================================================

z_0 = q_0 + 1i*p_0;
q = linspace(0,1,N); % q interval
p = linspace(-1/2,1/2,N); % p interval
[qmesh,pmesh]=meshgrid(q,p); 
z = (qmesh+1i*(pmesh)); 
Hus = zeros(N,N,t_final); % Average Husimi function array
cs=zeros(N,N);
norm_cs= (2/N)^0.25; % Normalisation constant for the coherent state
%==========================================================================
%   Create the initial state
%==========================================================================
% return

psi0=zeros(N,1);

for j = 1:N
   for m=-2:2
        psi0(j,1) = psi0(j,1) + exp(-pi*N*0.5*(abs(z_0)^2-(z_0)^2)-pi*N*(j/N-(z_0)+m)^2+(1i*pi*m)); 
   end
end


psi0=psi0*norm_cs;


% Projection onto the Schur basis

if ismember(Srep,'y')

psi0=psiS*psiS'*psi0;
psi0=psi0/norm(psi0);

end
% return
psi_2=zeros(N,t_final);

for t=1:t_final

    psi_2(:,t)=U^t*psi0;
    psi_2(:,t)=psi_2(:,t)./(psi_2(:,t)'*psi_2(:,t));
    
end
%==========================================================================
%   Create the split interval
%==========================================================================


n_split=2; % Who knew split has to be greater than 1
split=round(linspace(1,t_final,n_split));
psi_2_gpu=gpuArray(psi_2);
cs_gpu=gpuArray(cs);
q_gpu=gpuArray(q);
z_gpu=gpuArray(z);

% return
%==========================================================================
%   Begin split Husimi
%==========================================================================

for sn=2:length(split)
    
nstart=split(sn-1);
if sn>2 
    nstart=nstart+1; 
end
nend=split(sn);

% GPU Array Decleration
ds=abs(nend-nstart)+1;
Hus_split_gpu=gpuArray(zeros(N,N,ds));
psi_split_gpu=gpuArray(zeros(1,1,ds));

tic
for j = 1:N-1
    
   
    disp([num2str(j),' out of ',num2str(N),' for ',num2str(sn-1),' out of ',num2str(length(split)-1)]) % keep track
    cs_gpu=Cs_create_component(j,norm_cs,N,q_gpu,z_gpu,cs_gpu); 
    psi_split_gpu(1,1,:)=psi_2_gpu(j,nstart:nend);
    Hus_split_gpu=Hus_split_gpu+conj(cs_gpu).*psi_split_gpu;
    cs_gpu(:,:)=0;
    
end
toc

Hus(:,:,nstart:nend)=gather(Hus_split_gpu);
clear Hus_split_gpu psi_split_gpu

end
tic
% 
% n_efn=1

Hus_av=zeros(N,N);
for t=1:t_final
        disp([num2str(t),' out of ',num2str(t_final)]) % keep track 
        Hus_av=Hus_av+abs(Hus(:,:,t)).^2;
        
        figure(1)
        clf
        imagesc(q,p,abs(Hus(:,:,t)).^2)
        colorbar
        % caxis([0 1])
        colormap(viridis)
        set(gca,'YDir','normal')
        pause(0.5)

end
time2=toc;



viridis=viridis();

figure(2)
clf
imagesc(q,p,Hus_av./t_final)
colorbar
% caxis([0 1])
colormap(viridis)
set(gca,'YDir','normal')




% figure(5)
% clf
% imagesc(q,p,Hus_av)
% colormap(parula)
% set(gca,'YDir','normal')
% caxis([0 1])

% 
% if max(max(Hus_av))>1
%     'Hus_av has elements that are > 1!!!'
% end
% 
% 
% if sum(sum(isnan(Hus_av)))>0
%     'Hus_av has elements that are NaN!!!'
% end
% return
% while 1 
% 
% % Tell the user about the system and ask for input
%     
% user_msg_0 = strcat('Max Efns: ',num2str(n_efn));    
% display(user_msg_0);
% user_msg_1=' Select the number of eigenfunctions to look at ';
% nsplit=input(user_msg_1);
% 
% 
% Hus_av=zeros(N,N);
% for t=1:nsplit
%         disp([num2str(t),' out of ',num2str(n_efn)]) % keep track 
%         Hus_av=Hus_av+abs(Hus(:,:,t)).^2;
% end
% time2=toc;
% 
% 
% 
% viridis=viridis();
% 
% figure
% clf
% imagesc(q,p,Hus_av)
% colorbar
% % caxis([0 0.2])
% colormap(viridis)
% set(gca,'YDir','normal')
% 
% 
% user_exit=input(' Repeat with new nsplit? [y/n]:  ','s');
% 
% if isequal(user_exit,'n')
% 
%     break
%     
% end






% end



