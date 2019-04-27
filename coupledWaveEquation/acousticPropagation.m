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
n=5; m=n;
dx=.25; dy=dx;
fullDim=m*n;

% Populate the lower left quadrant $\mathbf{A_{21}}$, i.e. the discretized laplacian operator
r=[-2/dx^2-2/dy^2 1/dy^2 zeros(1, m-2) 1/dx^2 zeros(1,fullDim-(m+1))];
A21=toeplitz(r);

%%
% Set boundary element contributions to zero (artifact from making pressure
% matrix into a vector for state)
for i = 1:n-1 
    A21(i*n+1,i*n)=0;
    A21(i*n,i*n+1)=0;
end

%%
% Form the full A matrix
A=[zeros(n^2) eye(n^2);
   A21 zeros(n^2)];

%plot(svd(A))