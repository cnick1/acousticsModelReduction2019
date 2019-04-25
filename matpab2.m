%
% Build the matrices for the system.
%
function [a m iter]=matpab2(n,c,rho,dt,dx,dy)
p=ones(n^2-1,1);
for i=1:n-1
    p(i*n)=0;
end
a=[eye(n^2) zeros(n^2) dt/(rho*dx)*(eye(n^2)-diag(ones(n*(n-1),1),n));
   zeros(n^2) eye(n^2) dt/(rho*dy)*(eye(n^2)-diag(p,1));
   zeros(n^2) zeros(n^2) eye(n^2)];


vx=((dt*rho*c^2)/dx)*(eye(n^2)-diag(ones(n*(n-1),1),-n));
% vx([1:n],[1:n])=0;

vx([10.5*n-10:10.5*n+10],[1:end])=0;
vx([9.5*n-10:9.5*n+10],[1:end])=0;
vx([8.5*n-10:8.5*n+10],[1:end])=0;
vx([7.5*n-10:7.5*n+10],[1:end])=0;

vy=((dt*rho*c^2)/dy)*(eye(n^2)-diag(p,-1));
vy([10.5*n-10:10.5*n+10],[1:end])=0;
vy([9.5*n-10:9.5*n+10],[1:end])=0;
vy([8.5*n-10:8.5*n+10],[1:end])=0;
vy([7.5*n-10:7.5*n+10],[1:end])=0;

m=[eye(2*n^2) zeros(2*n^2,n^2); 
    vx vy eye(n^2)];
a=sparse(a);
m=sparse(m);
%m([2*n^2+1:2*n^2+n],[1:n])=zeros(n,n);
%m([2*n^2+n:2*n^2+n+5],[1:n])=0;
iter=m\a;