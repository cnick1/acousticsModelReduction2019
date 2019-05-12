%% Acoustic Wave Propagation using Wave Equation
% Model generation and reduction of acoustic wave propogation using the wave
% equation, i.e. 
% 
% $$ \ddot{p}=K \nabla^2 p $$
%
% where $p$ represents pressure and $K = \rho c^2$, where $\rho$ is the
% uniform density and $c$ is the speed of sound in the medium.
% In this version we will avoid state space and instead interpolate via
% Galerkin projection.
 
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
%% Formulating the system matrices $\mathbf{K}$ and $\mathbf{B}$
% Define grid size and spacing (square for now)
clear; close all; clc;
rho=1.2; c=345;
% rho=-1; c=1;

n=40; m=n;
dx=.25; dy=dx;
fullDim=m*n;

%%
% Populate the coefficient matrix $\mathbf{K}$, i.e. the discretized
% laplacian operator. It has an almost toeplitz form, with some extra zeros
% on the off diagonal for boundary element conditions.

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
% We also define our initial condition $\mathbf{B=p0}$, where
% $\mathbf{u}(t)$ is an impulse.

xx=1:1:n;
yy=1:1:n;
parab_x0=20;
parab_y0=parab_x0;
%choose parabHeight = height of parabola
parabHeight=.05;

parabRadius=10; % how many grid space for the radius of the parabola
beta=parabHeight/(parabRadius*dy)^2;
alfa=beta;
ff=zeros(n);
for i=1:n
    for j=1:n
        if -alfa*(xx(i) -parab_x0)^2-beta*(yy(j)-parab_y0)^2+parabHeight >= 0
        ff(i,j)=-alfa*(xx(i) -parab_x0)^2-beta*(yy(j)-parab_y0)^2+parabHeight; 
        end
    end
end
B=ff(:);

%% Solve the full model using ODE45 as a baseline
% Form the full A matrix and integrate using ODE45
% 
runFull=0;
if runFull
    A=[zeros(n^2) eye(n^2);
        -K zeros(n^2)];
    
    x0=zeros(1,2*m*n);
    x0(1:m*n)=B;
    
    tspan = [0 0.01];
    [t, p] = ode45(@(t,p) myfun(t,p,A), tspan, x0);
end

%% Form Galerkin reduction basis $\mathbf{V_r}$
% We form reduction basis $\mathbf{V_r}$ for the system, given by an
% orthonormal basis for the space spanned by the vectors 
%
% $\mathbf{V_i} = (s_i^2 \mathbf{I + K})^{-1} \mathbf{B}$
%

r=5^2;
si=logspace(1,5,r);

Vi=zeros(m*n,r);
I = eye(m*n);
for i=1:r
    Vi(:,i)=(si(i)^2*I+K)\B;
end

% figure
% loglog(si, abs(H(1,:)))

%%
% Rather than using Matlab's orth command, we use svd in case $\mathbf{V_i}$ is
% considered rank defficient. Orth calls svd according to Matlab
% documentation, so this just ensures that $\mathbf{V_i}$ has the correct dimensions
% and that $\mathbf{V_r^T V_r = I}$

[Vr,~] = svd(Vi,'econ');

%% Solve the reduced model using ODE45
% Set up the reduced system matrices using Galerkin projection
%
% $$ \mathbf{Kr = V_r^T K V_r} $$
%
% $$ \mathbf{Br = V_r^T B} $$
% 

Kr=Vr'*K*Vr;
Br=Vr'*B;

Ar=[zeros(r) eye(r);
   -Kr zeros(r)];

xr0=zeros(1,2*r);
xr0(1:r)=Br;

runReduced=0;
    if runReduced
    tspan = [0 .01];
    [t, pr] = ode45(@(t,pr) myfun(t,pr,Ar), tspan, xr0);
end

%% Plot the animations
%
if runFull==1
    fullP=p(:,1:m*n);

    for k = 1:size(fullP,1)
        surf(vectomat(fullP(k,:),m,n));
        axis([0 m 0 n -.1 .1])
        drawnow;
        Mframes(k) = getframe;
    end
end

if runReduced
    reducedP=Vr*pr(:,1:r)';

    for k = 1:size(reducedP, 2)
        surf(vectomat(reducedP(:,k),m,n));
        axis([0 m 0 n -.1 .1])
        drawnow;
        Mrframes(k) = getframe;
    end
end

function dp = myfun(t,p,A)
dp=A*p;
end