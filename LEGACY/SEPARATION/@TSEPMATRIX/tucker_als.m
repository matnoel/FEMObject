function [T,fit] = tucker_als(X,R,varargin)
%function [T,fit] = tucker_als(X,R,varargin)
%Higher-order orthogonal iteration.

[a,fit]=tucker_als(full(X.alpha),R,varargin{:});

a=TSEPMATRIX(a);
ga=gathervectors(a);
gx=gathervectors(X);

for i=1:X.dim
    ga.F{i}=gx.F{i}*ga.F{i};
end
T=splitvectors(ga);

end
