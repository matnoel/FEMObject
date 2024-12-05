function X = lhsrandom(rv,n,p)

if nargin==2
    oldn = n;
    n = p*n;
end

for i=1:length(rv.Y)
rv.Y{i} = lhsrandom(rv.Y{i},n,1);
end

for i=1:length(rv.funparam)
if isa(rv.funparam{i},'function_handle') | isa(rv.funparam{i},'inline')
    rv.funparam{i}=rv.funparam{i}(rv.Y{:});
end
end

X = rv.X(rv.funparam{:});
A = lhsdesign(n,1);
X = transfer(RVUNIFORM(0,1),X,A);

if nargin==2
    X =reshape(X,oldn,p);
end
