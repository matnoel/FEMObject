function c = ttv(a,v,dim)
%   c = ttv(a,v,dim)
%   tensor times vector along the dimension dim
%   Inspired by tensor/ttv from tensor toolbox

c = a;
% Permute it so that the dimensions we're working with come last
remdim = setdiff(1:ndims(a),dim);
if (ndims(a) > 1)
    c = permute(c,[remdim dim]);
end

n = ndims(a);
sz = size(a);
sz = sz([remdim dim]);

c=reshape(c,prod(sz(1:n-1)),sz(n));
c=c*v;
c=reshape(c,sz(1:n-1));
