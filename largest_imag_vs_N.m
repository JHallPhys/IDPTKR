clear all
close all


%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================

%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters
Nrange=[100,200,300,400,500,600,700,800,900,1000];
Krange=[1,2,2.5,5,10];
for ittk=1:length(Krange)
    ittk
    K_class = Krange(ittk);% Classical Kicking
    for ittn=1:length(Nrange)
        tic
        N_1 =Nrange(ittn); 
        N = 2*N_1+1; % Hilbert space dimension
        % return

        
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
        
        EIM(ittn,ittk)=max(EI);
        toc


    end
    
end


figure(1)
clf
hold on
plot(2*Nrange+1,EIM(:,1),'k.-','markersize',5)
plot(2*Nrange+1,EIM(:,2),'ko-','markersize',2)
plot(2*Nrange+1,EIM(:,3),'k+-','markersize',2)
plot(2*Nrange+1,EIM(:,4),'ksquare-','markersize',2)
plot(2*Nrange+1,EIM(:,5),'kdiamond-','markersize',2)
xlabel('N_1')
ylabel('$\displaystyle \tilde{\mu}_1$', 'Interpreter','latex')

figure(2)
clf
hold on
plot(Nrange,EIM(:,1),'k.-','markersize',5)
plot(Nrange,EIM(:,2),'ko-','markersize',5)
plot(Nrange,EIM(:,3),'k+-','markersize',5)
plot(Nrange,EIM(:,4),'ksquare-','markersize',5)
plot(Nrange,EIM(:,5),'kdiamond-','markersize',5)
xlabel('N')
ylabel('$\displaystyle \tilde{\mu}_1$', 'Interpreter','latex')
legend('k=1','k=2','k=2.5','k=5','k=10')
axis([min(Nrange) max(Nrange) 0 0.5])

