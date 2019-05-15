%
% Build the updating matrix A (U) for the system. a and m are the matrices
% with Px,Py and Vx,Vy respectively. See equation 5 in Tarazaga2008.
%
function [a m A]=matpab2(n,c,rho,dt,dx,dy)
p=ones(n^2-1,1);
for i=1:n-1
    p(i*n)=0;
end
a=[eye(n^2) zeros(n^2) dt/(rho*dx)*(eye(n^2)-diag(ones(n*(n-1),1),n));
   zeros(n^2) eye(n^2) dt/(rho*dy)*(eye(n^2)-diag(p,1));
   zeros(n^2) zeros(n^2) eye(n^2)];


Vx=((dt*rho*c^2)/dx)*(eye(n^2)-diag(ones(n*(n-1),1),-n));
Vy=((dt*rho*c^2)/dy)*(eye(n^2)-diag(p,-1));
%% Set any building conditions below

% Vx([1:n],[1:n])=0;
% Vx([10.5*n-10:10.5*n+10],[1:end])=0;
% Vx([9.5*n-10:9.5*n+10],[1:end])=0;
% Vx([8.5*n-10:8.5*n+10],[1:end])=0;
% Vx([7.5*n-10:7.5*n+10],[1:end])=0;
% 
% Vy([10.5*n-10:10.5*n+10],[1:end])=0;
% Vy([9.5*n-10:9.5*n+10],[1:end])=0;
% Vy([8.5*n-10:8.5*n+10],[1:end])=0;
% Vy([7.5*n-10:7.5*n+10],[1:end])=0;

%% Form full output matrices
m=[eye(2*n^2) zeros(2*n^2,n^2); 
    Vx Vy eye(n^2)];
a=sparse(a);
m=sparse(m);
A=m\a;