function Fi=margicdf(RF,X,x)
% function Fi=margicdf(RF,X,x)
% evaluate F_X^-1(x) : the marginal icdf of the random field RF
% X : coordinate of index  of size n * d (n number of evaluation points, d index dimension) 
% x : array
% Fi = ND-array of size  size(x,1) ... size(x,..) x n
pX=eval(RF.marginal,X);
keyboard
n=size(X,1);
m=size(x);
Fi=zeros(n,numel(x));
for i=1:numel(x)
    Fi(:,i)=icdf(pX,x(i));
end

Fi=reshape(Fi,[n,m]);
%Fi=squeeze(Fi);