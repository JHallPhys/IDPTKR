clear all
close all


%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================

%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters

N_1 =1000; 
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

% return
%==========================================================================
%   Project onto Subspace of stability
%==========================================================================

[psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values
Es=diag(Es);
[psi_2,n_efn]=Psi_lifetime(psiS,Es,eps,set_stability);
n_efn
return
ER=imag(log(diag(En)))/pi;
EI=real(log(diag(En)))/(imag(gamma)*N);

figure(1)
clf
plot(ER,EI,'k.','markersize',2)
% xlabel('\frac{\theta_n}{\pi}')
xlabel('$\displaystyle \tilde{\theta}_n$', 'Interpreter','latex')
ylabel('$\displaystyle \tilde{\mu}_n$', 'Interpreter','latex')
axis([0 1 -0.5 0.5])
return
%==========================================================================
%  Husimi Stuff 
%==========================================================================
psi_gamma=psiS(:,1:n_efn);
[q,p,z,dz]=get_husimi_grid(N);
Hus=get_husimi(N,n_efn,q,p,z,psi_gamma);
Hus_av=zeros(N,N);


for itt=1:n_efn
   
    Hus_av=Hus_av+abs(Hus(:,:,itt)).^2;
    
    if itt==1
    figure
    clf
    imagesc(q,p,abs(Hus(:,:,itt)).^2)
    % colorbar
    % caxis([0 1])
    colormap(viridis)
%     title(strcat('state',num2str(itt)))
    set(gca,'YDir','normal')
    end

end


str_title='Average_nefn_:';
str_title_num=num2str(n_efn);
str_gamma='_dgamma_';
str_gamma_num= strrep(num2str(imag(gamma)),'.','p');
str_average_title=strcat(str_title,str_title_num,str_gamma,str_gamma_num)

figure
clf
imagesc(q,p,Hus_av)
% colorbar
caxis([0 1])
colormap(viridis)
% title(str_average_title)
set(gca,'YDir','normal')

%==========================================================================
%  Entropy Stuff 
%==========================================================================

% return

Hus_Entropy=-Hus_av.*log(Hus_av);

for itt=1:n_efn
    Pj=abs(Hus(:,:,itt)).^2;
    Pj=Pj;
    Hus_Entropy=Pj;
    fname=fname_husimi_single_efn_special(K_class,N,gamma,itt,str_ext);
    parent_d = cd;  
    cd './Husimi_dat' % Directory where matrix is stored
    save(fname,'Hus_Entropy'); % save it 
    cd(parent_d)
end


















