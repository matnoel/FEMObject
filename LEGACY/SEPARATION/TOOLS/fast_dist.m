function d=fast_dist(a,b,varargin)
% Evalue rapidement la distance :
%    d(r) = || truncate(a,r) - b ||/||b||

if nargin == 2
    R=a.m;
elseif nargin ==3
    R=varargin{1};
end

d=zeros(1,min(a.m,R));
NB2=norm(b)^2;
dm1=1;
res=-b;
for r=1:min(a.m,R)
    a2=truncate(a,r);
    d(r)= dm1 + (2*fastprodscal( res,a2 ) + norm(a2)^2)/NB2;
    dm1=d(r);
    res=res+a2;
end

d=sqrt(abs(d));
