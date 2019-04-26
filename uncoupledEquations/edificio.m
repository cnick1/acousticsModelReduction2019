n=40;
xx=1:1:n;
yy=1:1:n;
x0=10;
x1=30
y0=7;
y1=10;
ff=zeros(n);
for i=x0:x1
    for j=y0:y1
        ff(i,j)=.001; 
    end
end
    
mesh(xx,yy,ff)
