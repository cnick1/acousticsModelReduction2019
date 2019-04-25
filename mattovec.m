%
%  This program write a vector out
%  of a matrix columnwise.
%      v=mattovec(x)
function v=mattovec(x)
[m,n]=size(x);
v=[];
for i=1:n
	v=[v;x(:,i)];
end
