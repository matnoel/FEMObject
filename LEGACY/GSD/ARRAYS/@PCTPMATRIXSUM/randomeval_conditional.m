function a = randomeval_conditional(x,i,xi,xj,varargin)

m = size(xj,2)+1;
noi = setdiff(1:m,i);
n = size(xj,1);
X = zeros(n,m);
if numel(xi)>1
error('') 
end
X(:,i)=xi;
X(:,noi)=xj;
xi = X;

a = randomeval(x,xi,varargin{:});
