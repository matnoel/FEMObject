function co = eval(correl,x1,x2)

n1 = size(x1);
if nargin==1
x2=zeros(1,n1(2));
end

n2 = size(x2);

L=correl.L;
if numel(L)==1
L = repmat(L,1,n1(2));
end

z=zeros(n1(1),n2(1));
for k=1:n1(2) 
x1k = repmat(x1(:,k),1,n2(1));
x2k = repmat(x2(:,k)',n1(1),1);
z=z+abs((x1k-x2k)/L(k));
end

co=exp(-z);