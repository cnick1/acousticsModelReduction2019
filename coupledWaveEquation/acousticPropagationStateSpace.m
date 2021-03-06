%% Acoustic Wave Propagation using Wave Equation
% Model generation and reduction of acoustic wave propogation using the wave
% equation, i.e. 
% 
% $$ \ddot{p}=K \nabla^2 p $$
%
% where $p$ represents pressure and $K = - \rho c^2$, where $\rho$ is the
% uniform density and $c$ is the speed of sound in the medium.
 
%% State-space representation via semi-discretization
% In order to perform model reduction, we first want to find the state
% space representation of our dynamical system, i.e. 
% 
% $$ \mathbf{\dot{x}}(t) = A \mathbf{x}(t) + B \mathbf{u}(t) $$
% 
% $$ \mathbf{y}(t) = C \mathbf{x}(t) + D \mathbf{u}(t) $$
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
% coupled ODEs, where the xy grid size is m by n. 
% 
% We can then define the state:
% 
% $$ \mathbf{x}(t) = \pmatrix{\mathbf{x}_1 \cr \mathbf{x}_2} = \pmatrix{\mathbf{p} \cr \mathbf{\dot{p}}}$$
% 
% where $\mathbf{p}$ is the columns of $p(i,j)$ (the pressure at each node)
% stacked to form a vector. Writing our equations out in state space form
% now, we obtain
%
% $$ \pmatrix{\mathbf{\dot{x}}_1 \cr \mathbf{\dot{x}}_2} = \pmatrix{\mathbf{0} & \mathbf{I} \cr \mathbf{A_{21}} & \mathbf{0}} \pmatrix{\mathbf{x}_1 \cr \mathbf{x}_2} $$
%
% where the equation for $\mathbf{\dot{x}_2}$ expands to the coupled system
% of ODEs for the pressure at each node. As such, $\mathbf{A_{21}}$ must
% multiply out to the finite difference formula for the laplacian in each
% equation.
%% Formulating the coefficient matrix $\mathbf{A}$
% Define grid size and spacing (square for now)
clear; close all; clc;
rho=1.2; c=345;
% rho=-1; c=1;

n=40; m=n;
dx=.25; dy=dx;
fullDim=m*n;

%%
% Populate the lower left quadrant $\mathbf{A_{21}}$, i.e. the discretized laplacian operator
r=[-2/dx^2-2/dy^2 1/dy^2 zeros(1, m-2) 1/dx^2 zeros(1,fullDim-(m+1))];
A21=rho * c^2 * toeplitz(r);

%%
% Set boundary element contributions to zero (artifact from making pressure
% matrix into a vector for state) (zero gives Dirichlet aka mirror BC).
for i = 1:n-1 
    A21(i*n+1,i*n)=0;
    A21(i*n,i*n+1)=0;
end

%%
% Form the full A matrix
A=[zeros(n^2) eye(n^2);
   A21 zeros(n^2)];

%plot(svd(A))
%% Solve the full model using ODE45
% Set up initial conditions of a parabola and integrate forward in time
% 

%%
% Initial condition setup

xx=1:1:n;
yy=1:1:n;
x0=20;
y0=20;
%choose C = height of parabola
C=.05;

k=15; % how many grid space for the radius of the parabola
beta=C/(k*dy)^2;
alfa=beta;
ff=zeros(n);
for i=1:n
    for j=1:n
        if -alfa*(xx(i) -x0)^2-beta*(yy(j)-y0)^2+C >= 0
        ff(i,j)=-alfa*(xx(i) -x0)^2-beta*(yy(j)-y0)^2+C; 
        end
    end
end
p0=zeros(1,2*m*n);
p0(1:m*n)=ff(:)';

tspan = [0 0.01];
[t, p] = ode45(@(t,p) myfun(t,p,A), tspan, p0);

for k = 1:size(p,1)
    lis(k)=sum(p(k,:));
end
for k = 1:size(p,1)
	surf(xx,yy,vectomat(p(k,:),m,n));
    axis([0 m 0 n -.1 .1])
    drawnow;
	Mframes(k) = getframe;
end



function dp = myfun(t,p,A)
dp=A*p;
end
