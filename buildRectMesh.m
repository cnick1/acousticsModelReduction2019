function [pos,vec]=buildRectMesh(numElements,sizeEle)

%function [pos,vec]=buildRectMesh(numElements,sizeEle)
%
%Builds a rectangular grid of elements with linked nodes
%
%[numElements] is a 2 element vector giving the number of elements in the x and then y directions
%[sizeEle] is a 2 element vector giving the size of the elements in the x and then y directions
%[pos] is a prod(numElements) by 2 matrix of x and y locations of the grid
%[vec] is a prod(numElements) by 4 matrix indication the adjacent elelments
%to the element in question the order is positive x element, negative x
%element, positive y element and negative y element. An elements number of
%prod(numElements)+1 implies a boundary
num_ele=prod(numElements);
x=0:sizeEle(1):(numElements(1)-1)*sizeEle(1);
y=0:sizeEle(2):(numElements(2)-1)*sizeEle(2);

X=repmat(x,1,numElements(2))';
Y=reshape(repmat(y,numElements(1),1),1,num_ele)';
pos=[X,Y];
% plot(pos(:,1),pos(:,2),'+')

%%
%calculating all of the adjacent nodes
vec=zeros(num_ele,4);
vec(:,1)=2:num_ele+1;
vec(:,2)=0:num_ele-1;
vec(:,3)=1+numElements(1):num_ele+numElements(1);
vec(:,4)=1-numElements(1):num_ele-numElements(1);
%fixing the edges
vec(numElements(1):numElements(1):num_ele,1)=(num_ele+1)*ones(numElements(2),1);
vec(1:numElements(1):num_ele,2)=(num_ele+1)*ones(numElements(2),1);
vec(num_ele-numElements(1)+1:num_ele,3)=(num_ele+1)*ones(numElements(1),1);
vec(1:numElements(1),4)=(num_ele+1)*ones(numElements(1),1);