%
%  This program write a matrix out
%  of a vector provided the number
%  of rows and columns.
%      x=vectomat(x,m,n)
function f=vectomat(x,m,n)
for i=1:n
	f(:,i)=x((i-1)*m+1:i*m);
end
f=f;