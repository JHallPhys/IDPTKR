clear all
close all


%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================

%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters
Nrange=[1250,1500,2000];
Krange=linspace(2*pi-pi/4,2*pi+pi/4,21)
% return
for itt_N=1:length(Nrange)
    itt_N
    tic
    for itt_kick=1:length(Krange)
        
        N_1 =Nrange(itt_N); 
        N = 2*N_1+1; % Hilbert space dimension
        % return
        K_class = Krange(itt_kick);% Classical Kicking
        
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
        
        EIM(itt_kick,itt_N)=max(EI);


    end
    toc
end


figure(1)
clf
hold on
plot(Krange,EIM(:,1),'r.-','markersize',2)
plot(Krange,EIM(:,2),'b.-','markersize',2)
plot(Krange,EIM(:,3),'g.-','markersize',2)
xlabel('k')
ylabel('$\displaystyle \tilde{\mu}_n$', 'Interpreter','latex')


