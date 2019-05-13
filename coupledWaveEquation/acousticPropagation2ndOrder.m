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
animate=1; %Set to 0 to save time when doing comments

rho=1.2; c=345;
% rho=-1; c=1;

xdim=40; ydim=xdim;
dx=.25; dy=dx;
n=xdim*ydim;
%%
% Populate the coefficient matrix $\mathbf{K}$, i.e. the discretized
% laplacian operator. It has an almost toeplitz form, with some extra zeros
% on the off diagonal for boundary element conditions.

row=[-2/dx^2-2/dy^2 1/dy^2 zeros(1, ydim-2) 1/dx^2 zeros(1,n-(ydim+1))];
K=-rho * c^2 * toeplitz(row);

%%
% Set boundary element contributions to zero (artifact from making pressure
% matrix into a vector for state) (zero gives Dirichlet aka mirror BC).
for i = 1:xdim-1 
    K(i*xdim+1,i*xdim)=0;
    K(i*xdim,i*xdim+1)=0;
end

%plot(svd(K))

%%
% We also define our initial condition $\mathbf{B=p0}$, where
% $\mathbf{u}(t)$ is an impulse.

xx=1:1:xdim;
yy=1:1:xdim;
parab_x0=20;
parab_y0=parab_x0;
%choose parabHeight = height of parabola
parabHeight=.05;

parabRadius=10; % how many grid space for the radius of the parabola
beta=parabHeight/(parabRadius*dy)^2;
alfa=beta;
ff=zeros(xdim);
for i=1:xdim
    for j=1:xdim
        if -alfa*(xx(i) -parab_x0)^2-beta*(yy(j)-parab_y0)^2+parabHeight >= 0
        ff(i,j)=-alfa*(xx(i) -parab_x0)^2-beta*(yy(j)-parab_y0)^2+parabHeight; 
        end
    end
end
B=ff(:);

%% Solve the full model using ODE45 as a baseline
% Form the full A matrix and integrate using ODE45
% 

A=[zeros(xdim^2) eye(xdim^2);
    -K zeros(xdim^2)];

x0=zeros(1,2*n);
x0(1:n)=B;

tf=0.01;
tspan = [0 tf];
[t, x] = ode45(@(t,p) myfun(t,p,A), tspan, x0);

fullP=x(:,1:n);


%% Form Galerkin reduction basis $\mathbf{V_r}$
% We form reduction basis $\mathbf{V_r}$ for the system, given by an
% orthonormal basis for the space spanned by the vectors 
%
% $\mathbf{V_i} = (s_i^2 \mathbf{I + K})^{-1} \mathbf{B}$
%

r=10;
irkaIter
si=logspace(1,5,r);

Vi=zeros(n,r);
I = eye(n);
for i=1:r
    Vi(:,i)=(si(i)^2*I+K)\B;
end

% figure
% loglog(si, abs(Vi(1,:)))

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

tspan = [0 tf];
[t, xr] = ode45(@(t,pr) myfun(t,pr,Ar), tspan, xr0);

% Project back to full space
reducedP=(Vr*xr(:,1:r)')';

% Interpolate reducedP to match the dimension of fullP
reducedPinterp=interp1(reducedP,linspace(1,size(reducedP,1),size(fullP,1)));

errorP=fullP-reducedPinterp;

%% Plot error for particular node
%
figure
plot(errorP(:,625))
hold on
plot(reducedPinterp(:,625))
plot(fullP(:,625))
%% Plot the animations
%

if animate  
    figure
    for k = 1:size(fullP,1)
        surf(vectomat(fullP(k,:)-reducedPinterp(k,:),ydim,xdim));
        axis([0 ydim 0 xdim -.1 .1])
        drawnow;
        errorFrames(k) = getframe;
    end
    
    for k = 1:size(fullP,1)
        surf(vectomat(fullP(k,:),ydim,xdim));
        axis([0 ydim 0 xdim -.1 .1])
        drawnow;
        Mframes(k) = getframe;
    end

    for k = 1:size(reducedPinterp, 1)
        surf(vectomat(reducedPinterp(k,:),ydim,xdim));
        axis([0 ydim 0 xdim -.1 .1])
        drawnow;
        Mrframes(k) = getframe;
    end
end

function dp = myfun(t,p,A)
dp=A*p;
end