%% Arnoldi model reduction on FDTD Acoustic Wave Propagation
% FDTD bench for model reduction
% 

%% Initialize constants and build updating matrix A
clear; close all; clc
rho=1.2; c=345;
n=40;
dt=.0001;
dx=.25;
dy=.25;
[~, ~, A]=matpab2(n, c, rho, dt, dx, dy); %Build the updating matrix A (U) for leapfrog technique
%% Initial condition setup
xx=1:1:n;
yy=1:1:n;
parabolaX0=20; % parabola center point
parabolaY0=20;
C=.05; % height of parabola

k=15; % how many grid space for the radius of the parabol
beta=C/(k*dy)^2;
alfa=beta;
ff=zeros(n);
for i=1:n
    for j=1:n
        if -alfa*(xx(i) -parabolaX0)^2-beta*(yy(j)-parabolaY0)^2+C >= 0
        ff(i,j)=-alfa*(xx(i) -parabolaX0)^2-beta*(yy(j)-parabolaY0)^2+C; 
        end
    end
end
p0=ff(:);
mesh(xx,yy,ff) % View initial condition

%%
ti=500;
% finaltime=dt*ti;
% lengthofsystem=n*dx;

fdtd0=[zeros(n^2,1);zeros(n^2,1);p0];   % for an impulse initial condition
%%
% Arnoldi
%break
redstep=100;
[V, H]=arnoldiMR(A,fdtd0,redstep);

% reduce A
Ah=H;
%reduce initial condition
fdtd0h=V'*fdtd0;

fdtdh=H*fdtd0h;
for i=2:ti
    fdtdh(:,i)=H*fdtdh(:,i-1);
end

fdtdhr=zeros(3*n*n,ti);
for i=1:ti
fdtdhr(:,i)=V*fdtdh(:,i);
end

cfdtdhr=fdtdhr(2*n^2+1:end,:); %Just the pressure variable from the reduced model
%Mhr(:,:,1)=vectomat(cfdtd(:,1),n,n);
% for i=1:ti
% Mhr(:,:,i)=vectomat(cfdtdhr(:,i),n,n);  % convert to matrix form
% end

%%
%fdtd=A*fdtd0+in*pul(1);
fdtd=A*fdtd0;   % impulse

for i=2:ti
    fdtd(:,i)=A*fdtd(:,i-1);  % iteration of the full size, real system
    %fdtd is a vector with all the states stacked up. They need to be
    %shaped into a matrix for a grid plot.
    
end

cfdtd=fdtd(2*n^2+1:end,:);

%%
x0=10;
x1=30;
y0=7;
y1=10;

M=zeros(n,n,ti);
for i=1:ti
M(:,:,i)=vectomat(cfdtd(:,i),n,n); % Convert pressure state vector to matrix for plotting

%     for ii=x0:x1 %This shows the building
%         for jj=y0:y1
%         M(ii,jj,i)=.01;
%         end
%     end
end



%movie(F,1,2)
jj=1;
for i=1:5:ti

surf((1:n),(1:n),M(:,:,i))

view(-28, 66)
%axis([0 10 0 10 -.05 .05])
F(jj)=getframe(1,[38,30,473,373]);
jj=jj+1;
end
%%
% reduced case

Mhr=zeros(n,n,ti);
for i=1:ti
Mhr(:,:,i)=vectomat(cfdtdhr(:,i),n,n); % convert to matrix form
%     for ii=x0:x1 %This shows the building
%         for jj=y0:y1
%         Mhr(ii,jj,i)=.01; 
%         end
%     end
%fmodel(i)=M(15,20,i);
end



%movie(F,1,2)
jj=1;
for i=1:6:ti

surf((1:n),(1:n),Mhr(:,:,i))

view(-28, 66)
%axis([0 10 0 10 -.05 .05])
G(jj)=getframe(1,[38,30,473,373]);
jj=jj+1;
end
%%
figure
% Time plot comparisons
point=11.5*n;
rmodel=cfdtdhr(point,:);
fmodel=cfdtd(point,:);
plot(1:ti,fmodel,'b')
axis([0 500 -10e-3 13e-3])
hold on
plot( 1:ti,rmodel,'r:')
xlabel('time steps')
ylabel('Pressure')
grid on
legend( 'FDTD',['r=',num2str(redstep),', Arnoldi'])
