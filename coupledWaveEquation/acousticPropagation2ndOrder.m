%% Acoustic Wave Propagation using Wave Equation
% Model generation and reduction of acoustic wave propogation using the wave
% equation, i.e. 
% 
% $$ \ddot{p}=K \nabla^2 p $$
%
% where $p$ represents pressure and $K = - \rho c^2$, where $\rho$ is the
% uniform density and $c$ is the speed of sound in the medium.
% In this version we will avoid state space and instead interpolate the
% transfer function.
 
%% Semi-discretization
% In order to perform model reduction, we first want to develop a model for
% our system which looks like:
% 
% $$ \mathbf{M \ddot{p}} + \mathbf{ D \dot{x}} + \mathbf{ K p} = \mathbf{B u}(t) $$
%  
% We start by semi-discretizing the wave equation in space; this transforms
% the laplacian, which contains spacial first derivatives, into a
% collection of finite differences, i.e. 
% 
% $$ \nabla^2 p(i,j) = \frac{p(i+1,j)-2 p(i,j)+p(i-1,j)}{(\Delta x)^2} +
% \frac{p(i,j+1)-2 p(i,j)+p(i,j-1)}{(\Delta y)^2} $$
% 
% <matlab:web('https://en.wikipedia.org/wiki/Finite_difference_method#Example:_The_Laplace_operator','-browser') Laplacian Finite Difference Wikipedia>
% 
% This transforms the wave equation from a 2nd order PDE into m*n 2nd order
% coupled ODEs, where the xy grid size is m by n. The laplacian is turned
% into the coefficient matrix $\mathbf{K}$, which can be deduced from the
% discretization shown above.
% 
%% Formulating the coefficient matrix $\mathbf{K}$
% Define grid size and spacing (square for now)
clear; close all; clc;
rho=1.2; c=345;
% rho=-1; c=1;

n=40; m=n;
dx=.25; dy=dx;
fullDim=m*n;

%%
% Populate the coefficient matrix $\mathbf{K}$, i.e. the discretized laplacian operator
row=[-2/dx^2-2/dy^2 1/dy^2 zeros(1, m-2) 1/dx^2 zeros(1,fullDim-(m+1))];
K=-rho * c^2 * toeplitz(row);

%%
% Set boundary element contributions to zero (artifact from making pressure
% matrix into a vector for state) (zero gives Dirichlet aka mirror BC).
for i = 1:n-1 
    K(i*n+1,i*n)=0;
    K(i*n,i*n+1)=0;
end

%plot(svd(K))

C=eye(m*n);

%%
% We also define our initial condition $\mathbf{p0}$

xx=1:1:n;
yy=1:1:n;
x0=10;
y0=x0;
%choose parabHeight = height of parabola
parabHeight=.05;

k=2; % how many grid space for the radius of the parabola
beta=parabHeight/(k*dy)^2;
alfa=beta;
ff=zeros(n);
for i=1:n
    for j=1:n
        if -alfa*(xx(i) -x0)^2-beta*(yy(j)-y0)^2+parabHeight >= 0
        ff(i,j)=-alfa*(xx(i) -x0)^2-beta*(yy(j)-y0)^2+parabHeight; 
        end
    end
end
p0=ff(:);

%% Form Transfer Function $\mathbf{H}(s)$
% We form and plot the transfer function for the system, given by 
%
% $\mathbf{H}(s) = \mathbf{C}(s^2 + \mathbf{K})^-1 \mathbf{p0}$
%

sampleNum=16^2;
omegas=logspace(1,5,sampleNum);

H=zeros(m*n,sampleNum);
I = eye(m*n);
for i=1:sampleNum
    H(:,i)=((1i*omegas(i))^2*I+K)\p0;
end



figure
loglog(omegas, abs(H(1,:)))

r=sampleNum;
si=omegas;
Vr=H;

Kr=Vr'*K*Vr;



%% Solve the full model using ODE45
% Set up initial conditions of a parabola and integrate forward in time
% 

%%
% Form the full A matrix
A=[zeros(n^2) eye(n^2);
   -K zeros(n^2)];

pr0=zeros(1,2*m*n);
pr0(1:m*n)=ff(:);

tspan = [0 0.01];
[t, p] = ode45(@(t,p) myfun(t,p,A), tspan, pr0);

for k = 1:size(p,1)
    lis(k)=sum(p(k,:));
end
for k = 1:size(p,1)
	surf(vectomat(p(k,:),m,n));
    axis([0 m 0 n -.1 .1])
    drawnow;
	Mframes(k) = getframe;
end


%% Solve the reduced model using ODE45
% Set up initial conditions of a parabola and integrate forward in time
% 

Ar=[zeros(r) eye(r);
   -Kr zeros(r)];

pr0=zeros(1,2*r);
pr0(1:r)=Vr'*p0;

tspan = [0 .01];
[t, pr] = ode45(@(t,pr) myfun(t,pr,Ar), tspan, pr0);

for k = 1:size(pr,1)
    lis(k)=sum(pr(k,:));
end
for k = 1:size(pr,1)
	surf(vectomat(pr(k,:),sqrt(r),sqrt(r)));
    axis([0 sqrt(r) 0 sqrt(r) -.1 .1])
    drawnow;
	Mrframes(k) = getframe;
end



function dp = myfun(t,p,A)
dp=A*p;
end