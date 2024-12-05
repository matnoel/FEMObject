function r = randomeval(r,xi,varargin)
% function r = randomeval(r,xi,varargin)

if isa(r,'cell')
    for l=1:length(r)
        r{l} = randomeval(r{l},xi,varargin{:});
    end
elseif isa(r,'FunctionalTensor') || isa(r,'FunctionalBasisArray')
    r = r.eval(xi,varargin{:});
else
    if size(xi,1)==1
        r = r;
    elseif numel(r)==1
        r = repmat(r,1,size(xi,1));
    else
        s = size(r);
        sm = [size(xi,1),1];
        r = repmat(r(:),1,size(xi,1));
        r = MULTIMATRIX(r,s,sm);
    end
end