clear all;
close all

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Hi Eva :) 
% As a reminder
% The Interesting Quantum Stuff happens for gamma=0.001, gamma=0.005
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%==========================================================================


%==========================================================================
%   System Parameters 
%==========================================================================

k=2*pi
% k=1.1
gamma =0.001;  % PT Strength
tmax =100; % % Number of itterations
N = 25; % Number of Initial conditions (N^2)

%==========================================================================
%   Grid of Initial Conditions
%==========================================================================

% dN=1/N^2


% qmin=1;
% qmax=0;    
% pmin=-0.5;
% pmax=0.5;

qmin=0.2325;
qmax=0.255;    
pmin=-6*gamma;
pmax=6*gamma;
% 


% 
q=linspace(qmin,qmax,N);
p = linspace(pmin,pmax,N);
% 

% q=linspace(0,1,N);
% p = linspace(-0.5,0.5,N);


[qmesh,pmesh] = meshgrid(q,p);
[qmesh_0,pmesh_0] = meshgrid(q,p);




%==========================================================================
%   Dynamics
%==========================================================================
norm_p=zeros(tmax,1);
norm_p0=1;
for t=1:tmax % itteration time

t
    
% qmesh_old=qmesh;
% qmesh=mod(qmesh+pmesh+k*sin(2*pi*qmesh)/(4*pi),1);
% pmesh=mod(pmesh+k*(sin(2*pi*qmesh_old)+sin(2*pi*qmesh))/(4*pi)+1/2,1)-1/2;    
%     

figure(1)
clf
hold on
plot(qmesh,pmesh,'k.','Markersize',1)
xlabel('q')
ylabel('p')
axis([0 1 -0.5 0.5])
set(gca,'FontSize',10)
title('Start')
% pause(0.1)

% pstar=pmesh+k*sin(2*pi*qmesh)/(4*pi);
% qmesh=mod(qmesh+pstar+gamma/(2*pi),1);
% pmesh=mod(pstar+k*sin(2*pi*qmesh)/(4*pi)+gamma/(2*pi)+1/2,1)-1/2;    
% figure(2)
% clf
% plot(qmesh,pmesh,'k.','Markersize',1)
% xlabel('q')
% ylabel('p')
% axis([0 1 -0.5 0.5])
% set(gca,'FontSize',10)
% title('FWD')
% pstar=pmesh+k*sin(2*pi*qmesh)/(4*pi);
% qmesh=mod(qmesh+pmesh+k*sin(2*pi*qmesh)/(4*pi)+gamma/(2*pi),1);
% pmesh=mod(pstar+k*sin(2*pi*qmesh)/(4*pi)+1/2,1)-1/2;    
% 
% figure(1)
% clf
% plot(qmesh,pmesh,'k.','Markersize',1)
% xlabel('q')
% ylabel('p')
% axis([0 1 -0.5 0.5])
% set(gca,'FontSize',10)
% pause(2)

pstar=pmesh-k*sin(2*pi*qmesh)/(4*pi);
qmesh=mod(qmesh-pstar+gamma/(2*pi),1);
pmesh=mod(pstar-k*sin(2*pi*qmesh)/(4*pi)+gamma/(2*pi)+1/2,1)-1/2;    
   
% figure(3)
% clf 
% plot(qmesh,pmesh,'k.','Markersize',1)
% xlabel('q')
% ylabel('p')
% axis([0 1 -0.5 0.5])
% set(gca,'FontSize',10)
% pause(1)
% title('BWD')
% return
end


% figure(2)
% hold on
% plot(1:1:t,p_av,'k.-','Markersize',5)

% figure(3)
% hold on
% plot(1:1:t,norm_p,'k.-','Markersize',5)

% 
% figure(1)
% hold on 
% plot(qmesh,pmesh,'k.','Markersize',1)
