clear ; clc; close all;
%Physical parameters
mu=1.8*1e-3;
rou=0.917*1e3;
lf=3.34*1e5;
gama=3.2*1e-2;
Tm=273.15;
V=7.33*1e-6;
G=9/4*1e2;
np=0.16;
a=np*1e-6;
tau=4*gama*Tm/a/rou/lf/G/V;
A=3.65*1e-21;
phi0=0.2514;
phip=0.637;
kepa=a^2*(1-phip)^3/45/phip^2;
phi=phi0/(phip-phi0);
M=mu*V*Tm/kepa/G/rou/lf;
Mphi=M*(phip-phi0)/phip;
N=48*3.14*mu*a^2*V/A;
nt=2;
v0=2/(sqrt(1+4*M*phi)+1);
ttran=(M*phi*v0/N/(1-v0)^3)^(1/2);
%ttran2=((M*phi*v0+v0-1)/N/(1-v0)^3)^(1/2);
dt=0.002;
t=0:dt:nt;
h=zeros(nt/dt+1,1);
z=zeros(nt/dt+1,1);
%vi=zeros(nt/dt+1,1);

% Setup initial position profile
vi(1,1)=v0;
h(1,1)=0;
z(1,1)=0;
h(2,1)=vi(1,1)*phi*dt;
z(2,1)=(vi(1,1)-1)*dt;
% Set boundary conditions
cnt=1;
for i=2:nt/dt % Timestep loop
% Compute new position
z(i+1,1)=(z(i,1)/dt-1)/(1/dt+1/(M*h(i,1)-N*(z(i,1))^3));
h(i+1,1)=h(i,1)+dt*phi*(-z(i+1,1)/(M*h(i,1)-N*(z(i+1,1))^3));
%vi(i+1,1)=-z(i+1,1)/(M*h(i+1,1)-N*(z(i+1,1))^3);
%h(i+1,1)=h(i,1)+dt*phi*vi(i+1,1);
%z(i+1,1)=z(i,1)+dt*(-1+vi(i+1,1));

if(mod(i,50) == 0)
    % Plot solution
    figure(1), clf
    plot(t(1:i+1),h(1:i+1,1),'b');
    hold on;
    plot(t(1:i+1),-z(1:i+1,1),'r');
    xlabel('time (time scale)')
    ylabel('position (length scale)')
    drawnow
    %errorbar(tr1,pr,e,'or','markersize',6,'markerfacecolor','red')
end



if(abs(z(i,1)+h(i,1))<0.0003) & z(i) ~ 0
    data1(cnt) = i;
    data2(cnt,:) = [-z(i),h(i)]; cnt =cnt + 1;
end
end