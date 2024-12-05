function w = concat(u,v,k)
% function w = concat(u,v,k)

nu = [size2D(u),sizeND(u)];
nv = [size2D(v),sizeND(v)];

nu = [nu,ones(1,max(length(nv)-length(nu),0))];
nv = [nv,ones(1,max(length(nu)-length(nv),0))];

if nargin<3
    k=3;
end

if k>length(nu)
    nu = [nu,ones(1,k-length(nu))];
end
if k>length(nv)
    nv = [nv,ones(1,k-length(nv))];
end

nku = setdiff(1:ndims(nu),k);
nkv = setdiff(1:ndims(nv),k);

if any(nu(nku)~=nv(nkv))
    error('pour etre concatener suivant une dimension deux matrices, les autres dimensions doivent correspondre')
end

nw = nu;
nw(k) = nu(k)+nv(k);
w = u;
w.double = zeros(nw);
rep1 = cell(1,length(nw));
rep1(:) = {':'};
rep1{k} = 1:nu(k);
w.double(rep1{:}) = u.double;
rep1{k} = nu(k)+1:nu(k)+nv(k);
w.double(rep1{:}) = v.double;
