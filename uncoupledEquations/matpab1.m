%
%  matrix of the system
%
function a =matpab1(n,c,rho,dt,dx,dy)
p=ones(n^2-1,1);
for i=1:n-1
    p(i*n)=0;
end
a=[eye(n^2) zeros(n^2) (dt/(rho*dx)*(eye(n^2)-diag(ones(n*(n-1),1),n));
   zeros(n^2) eye(n^2) (dt/(rho*dy)*(eye(n^2)-diag(p,1));
   ((dt*rho*c^2)/dx)*(-eye(n^2)+diag(ones(n*(n-1),1),-n)) ((dt*rho*c^2)/dy)*(-eye(n^2)+diag(p,-1)) eye(n^2)];

