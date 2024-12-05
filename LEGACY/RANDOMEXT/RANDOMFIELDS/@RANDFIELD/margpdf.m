function p=margpdf(RF,X,x)
% function p=margpdf(RF,X,x)
% evaluate p_X(x) : the marginal pdf of the random field RF
% X : coordinate of index  of size n * d (n number of evaluation points, d index dimension) 
% x : array
% p = ND-array of size  size(x,1) ... size(x,..) x n
pX=eval(RF.marginal,X);
n=size(X,1);
m=size(x);

p=zeros(n,numel(x));
for i=1:numel(x)
    p(:,i)=pdf(pX,x(i));
end

p=reshape(p,[n,m]);
%p=squeeze(p);