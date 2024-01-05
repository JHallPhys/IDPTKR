%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================
clear all
close all
%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters
% Nrange=[250,500,750,1000];
% Nrange=[1250,1500,1750,2000];
Nrange=[250,500,750,1000,1250,1500,1750,2000];
% Nrange=[250,500,750];
Nbinz=2*max(Nrange);
Nedges=linspace(-0.5,0.5,Nbinz);
Pstate=zeros(Nbinz-1,length(Nrange));
Pstate_normalised=zeros(Nbinz-1,length(Nrange));
for ittn=1:length(Nrange)
ittn
tic

N_1 =Nrange(ittn); 
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

ER=imag(log(diag(En)))/pi;
EI=real(log(diag(En)))/(imag(gamma)*N);


%==========================================================================
%  Make the distribution of eigenvalues 
%==========================================================================

if ittn==1 % Get the maximum for the figure 
    frac_max=1-2*n_efn/N;
end


for ittstate=1:Nbinz-1
%     find(EI>Nedges(ittstate) && EI<Nedges(ittstate+1) )
%     find(EI>Nedges(ittstate) && EI<Nedges(ittstate+1) )
%     Pstate(ittstate,ittn)=length(find(EI>Nedges(ittstate) & EI<=Nedges(ittstate+1) ));
        Pstate(ittstate,ittn)=length(find(EI>Nedges(ittstate)))/N;
        Pstate_normalised(ittstate,ittn)=length(find(EI./max(EI)>2*Nedges(ittstate)))/N; % Plot the normalised distribution
end
% Pstate(:,ittn)=Pstate(:,ittn)/(N);
toc
end

Gain_range=linspace(-0.5,0.5,Nbinz-1);
Gain_range=Gain_range(round(Nbinz/2)+1:Nbinz-1);

figure(1)
clf
hold on
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,1),'r.-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,2),'b.-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,3),'g.-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,4),'m.-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,1),'rdiamond-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,2),'b.-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,3),'g.-','markersize',1,'Linewidth',1)
plot(linspace(-0.5,0.5,Nbinz-1),Pstate(:,4),'m.-','markersize',1,'Linewidth',1)
% plot(linspace(-0.5,0.5,Nbinz-1),Pstate./N,'g.-','markersize',2)
xlabel('$\displaystyle \mu$', 'Interpreter','latex')
ylabel('$\displaystyle P(\tilde{\mu})$', 'Interpreter','latex')

figure(2)
clf
hold on
plot(Gain_range,Pstate(round(Nbinz/2)+1:Nbinz-1,1),'r.-','markersize',1,'Linewidth',1)
plot(Gain_range,Pstate(round(Nbinz/2)+1:Nbinz-1,2),'b.-','markersize',1,'Linewidth',1)
plot(Gain_range,Pstate(round(Nbinz/2)+1:Nbinz-1,3),'g.-','markersize',1,'Linewidth',1)
plot(Gain_range,Pstate(round(Nbinz/2)+1:Nbinz-1,4),'m.-','markersize',1,'Linewidth',1)
% plot(linspace(-0.5,0.5,Nbinz-1),Pstate./N,'g.-','markersize',2)
xlabel('$\displaystyle \mu$', 'Interpreter','latex')
ylabel('$\displaystyle P(\tilde{\mu})$', 'Interpreter','latex')
% legend('N_1=250','N_1=500','N_1=750','N_1=1000')

% 
% figure(3)
% clf
% hold on
% plot(2*linspace(-0.5,0.5,Nbinz-1),Pstate_normalised(:,1),'r.-','markersize',1,'Linewidth',1)
% plot(2*linspace(-0.5,0.5,Nbinz-1),Pstate_normalised(:,2),'b.-','markersize',1,'Linewidth',1)
% plot(2*linspace(-0.5,0.5,Nbinz-1),Pstate_normalised(:,3),'g.-','markersize',2,'Linewidth',1)
% plot(2*linspace(-0.5,0.5,Nbinz-1),Pstate_normalised(:,4),'m.-','markersize',2,'Linewidth',1)
% % plot(linspace(-0.5,0.5,Nbinz-1),Pstate./N,'g.-','markersize',2)
% xlabel('$\displaystyle \mu$', 'Interpreter','latex')
% ylabel('$\displaystyle P({\mu})$', 'Interpreter','latex')

% figure(3)
% clf
% hold on
% plot(linspace(-0.5,0.5,Nbinz-1),0.5-abs(0.5-Pstate(:,1)),'r.-','markersize',2,'Linewidth',1)
% plot(linspace(-0.5,0.5,Nbinz-1),0.5-abs(0.5-Pstate(:,2)),'b.-','markersize',2,'Linewidth',1)
% plot(linspace(-0.5,0.5,Nbinz-1),0.5-abs(0.5-Pstate(:,3)),'g.-','markersize',2,'Linewidth',1)
% plot(linspace(-0.5,0.5,Nbinz-1),0.5-abs(0.5-Pstate(:,4)),'m.-','markersize',2,'Linewidth',1)
% % plot(linspace(-0.5,0.5,Nbinz-1),Pstate./N,'g.-','markersize',2)
% xlabel('$\displaystyle (\mu)$', 'Interpreter','latex')
% ylabel('$\displaystyle \tilde{P}({\mu})$', 'Interpreter','latex')
% axis([-0.5 0.5 0 frac_max])
% return