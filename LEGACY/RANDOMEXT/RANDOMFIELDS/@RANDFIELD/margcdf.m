function F=margcdf(RF,X,x)
% function F=margcdf(RF,X,x)
% evaluate F_X(x) : the marginal cdf of the random field RF
% X : coordinate of index  of size n * d (n number of evaluation points, d index dimension) 
% x : array
% F = ND-array of size  size(x,1) ... size(x,..) x n
pX=eval(RF.marginal,X);
n=size(X,1);
m=size(x);
F=zeros(n,numel(x));
for i=1:numel(x)
    F(:,i)=cdf(pX,x(i));
end

F=reshape(F,[n,m]);
F=squeeze(F);