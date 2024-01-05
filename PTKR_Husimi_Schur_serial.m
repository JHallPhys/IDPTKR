clear all
close all


%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================

%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters

N_1 =100; 
N = 2*N_1+1; % Hilbert space dimension
% return
K_class = 10% Classical Kicking

% cool at k=12.716,14.125
% K_class = 7.54545 % Classical Kicking
% K_class=12.6;
gamma = complex(0,0.001); % PT-strength 
T=2*pi/N; % Effective hbar
kick = K_class/T; % Quantum Kicking
hbar_eff=1/(2*pi*N);
str_ext='.mat';
% Schur Parameters

% eps=exp(imag(gamma)); % Tolerance parameter for stability classification
eps=1+1e-6;
% eps=1+hbar_eff/2
set_efn='G'; % Invarient Subspace: Gain ('G'), Stable ('S'), Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Stable ('0'), Loss ('-')

% Check that valid set_efn and set_stability labels are used

% return
if ismember(set_efn,{'G','L'})~=1 % Sorting label check
    display('Bad Eigenfunction set label')
    return
end


if ismember(set_stability,{'+','-','0'})~=1 % Subspace label check
    display('Bad set stability label')
    return
end


% return
%==========================================================================
%   Matrix Construction and Schur decompesition
%==========================================================================

U=zeros(N,N); % Initialise Flouqet matrix
U=UCheck(N,N_1,K_class,T,gamma,str_ext);% Check if matrix exists, if it does load it, else make and save it
[psi,En]=ECheck(U,N,N_1,K_class,T,gamma,str_ext);% Check if matrix exists, if it does load it, else make and save it
%==========================================================================
%   Project onto Subspace of stability
%==========================================================================

[psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values
[psi_2,n_efn]=Psi_lifetime(psiS,Es,eps,set_stability);
% 
% for j=1:n_efn
%     
%     G=real(log(psi_2(:,j)'*U*psi_2(:,j)));
%     
%     figure(2)
%     hold on 
%     plot(j,G,'b.','markersize',1)
% 
% end
% 
% 
% return



n_efn
return
% % return
n_efn/N
% return


% ylim([-0.15 0.15])
% ylabel('Im(\epsilon)')
% xlabel('Re(\epsilon)/ \pi')
% return 
%==========================================================================
%  Arrays independent of splitying
%==========================================================================

% Discrete Phase Space Grid
n_efn=1
q = linspace(0,1,N); % q interval
p = linspace(-1/2,1/2,N); % p interval
[qmesh,pmesh]=meshgrid(q,p); 
z = (qmesh+1i*(pmesh)); 
Hus = zeros(N,N); % Average Husimi function array
cs=zeros(N,N);
norm_cs= (2/N)^0.25; % Normalisation constant for the coherent state
efn_pick=1;
psi=psi_2(:,efn_pick);
tic
for j = 1:N-1
    j
   
    cs=Cs_create_component(j,norm_cs,N,q,z,cs); 
    Hus=Hus+conj(cs).*psi(j);
    cs(:,:)=0;
    
end
toc

Hus_av=abs(Hus).^2;


viridis=viridis();

figure
clf
imagesc(q,p,Hus_av)
colorbar
% caxis([0 1])
colormap(viridis)
set(gca,'YDir','normal')

Hus_Entropy=-Hus_av.*log(Hus_av);

figure
clf
imagesc(q,p,Hus_Entropy)
colorbar
% caxis([0 1])
colormap(viridis)
set(gca,'YDir','normal')

fname=fname_husimi_single_efn_special(K_class,N,imag(gamma),1,str_ext);
parent_d = cd;  
cd './Husimi_dat' % Directory where matrix is stored
save(fname,'Hus_Entropy'); % save it 
cd(parent_d)



if max(max(Hus_av))>1
    'Hus_av has elements that are > 1!!!'
end


if sum(sum(isnan(Hus_av)))>0
    'Hus_av has elements that are NaN!!!'
end










