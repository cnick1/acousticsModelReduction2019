n=40;
xx=1:1:n;
yy=1:1:n;
x0=10;
y0=10;
%chosse C = height of parabol
C=1;
dy=1;
k=3; % how many grid space for the radius of the parabol
beta=C/(k*dy)^2;
alfa=beta;
ff=zeros(n);
for i=1:n
    for j=1:n
        if -alfa*(xx(i) -x0)^2-beta*(yy(j)-y0)^2+C >= 0;
        ff(i,j)=-alfa*(xx(i) -x0)^2-beta*(yy(j)-y0)^2+C; 
        end
    end
end
    
%surf(xx,yy,ff)
%figure
mesh(xx,yy,ff)
