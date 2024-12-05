function a = random(apc,varargin)
% function a = random(apc,varargin)

[rep,pos] = isclassin('double',varargin);
rep = ischarin('init',varargin);
if rep
    initstate
end
n = prod(varargin{pos});
D = random(apc.POLYCHAOS,n);
switch ndims(apc)
    case 2
        a = double(apc.MYDOUBLE);
        if apc.pcdim==1
            a = D*a;
        else
            a = a*D';
        end
    otherwise
        nopc = setdiff(1:ndims(apc),apc.pcdim);
        a = double(apc.MYDOUBLE);
        a = permute(a,[nopc , apc.pcdim]);
        siz = size(a);
        a = reshape(a,prod(siz(1:end-1)),siz(end));
        a = a*D';
        a = reshape(a,[siz(1:end-1),n]);
        rep = 1:ndims(a);
        rep(nopc) = 1:ndims(a)-1;
        rep(apc.pcdim) = ndims(a);
        a = permute(a,rep);
end