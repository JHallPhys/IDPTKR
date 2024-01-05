clear all
close all
%==========================================================================
%   Initialisation  
%==========================================================================
set_efn='G'; % % Thou shalt not change
set_stability='+'; % Thou shalt not change
K_class =1.1; % Classical Kicking
gamma = complex(0,0.001); % PT-strength 
Nmin=200;
dN=200;
Nmax=1000;
N_i=Nmin:dN:Nmax;
NLen=length(N_i);
str_ext='.mat'; % File extension
% return
%==========================================================================
%  Fractal Weyl
%==========================================================================
for itt=1:NLen
    tic
    N_1=N_i(itt)
    N = 2*N_1+1; % Hilbert space dimension
    T=2*pi/N; % Effective hbar
    kick = K_class/T; % Quantum Kicking
    hbar=1/(N);
    % Schur Parameter
    eps=exp(imag(gamma)*T);% Tolerance parameter for stability classification
    U=UCheck(N,N_1,K_class,T,gamma,str_ext);% Check if matrix exists, if it does load it, else make and save it
    [psi,En]=ECheck(U,N,N_1,K_class,T,gamma,str_ext);% Check if matrix exists, if it does load it, else make and save it
    [psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values
    Es=diag(Es); 
    [psi_2,n_efn]=Psi_lifetime(psiS,Es,eps,set_stability);
    Es=1i*log(Es);
    pint(itt)=n_efn;
    toc
%     return 
end
%==========================================================================
% save the data  
%==========================================================================

FWData=[transpose(2*N_i+1),transpose(pint)]; % store the data
% Savce the data
% fname1=get_fname_ptkr(Nmin,dN,Nmax,K_class,gamma,set_stability,str_ext)
% parent_d = cd;   
% cd '/Users/joe/Desktop/PTFW/PTKR_FWeyl/FWdata'; % Where save data
% save(fname1,'FWData');
% cd(parent_d)


coeff = polyfit(log10(2*N_i+1),log10(pint), 1);
dfit = polyval(coeff , log10(2*N_i+1));


figure
hold on
plot(log10(2*N_i+1),log10(pint),'k.-','Markersize',10)
plot(log10(2*N_i+1),dfit,'r.-','Markersize',10)
xlabel('N')
ylabel('pint(\mu)')
FWDimension=coeff(1)