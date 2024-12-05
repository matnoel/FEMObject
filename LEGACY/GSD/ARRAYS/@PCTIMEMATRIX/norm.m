function an = norm(a,varargin)
% function an = norm(a,varargin)

if isa(a.value,'PCRADIALMATRIX')
    an = norm(a.value,[],getMmatrix(a));
else
    
    an = full(abs(prodscal(a,a,varargin{:})));
    
    if an<100*eps
        an = an;
    else
        an = sqrt(an);
    end
end
