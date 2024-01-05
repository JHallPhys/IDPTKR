function [q,p,z,dz]=get_husimi_grid(N)

q = linspace(0,1,N); % q interval
p = linspace(-0.5,0.5,N); % p interval
dq=abs(q(2)-q(1));
dp=abs(p(2)-p(1));
[qmesh,pmesh]=meshgrid(q,p); 
z = (qmesh+1i*(pmesh)); 
dz=dq*dp;
end