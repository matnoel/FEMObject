function am=mean(apc)

Mfe = findpolyfe(apc);

if isempty(Mfe)
n=ndims(apc.MYDOUBLE);
rep = cell(1,n);
rep(:)={':'};
rep(apc.pcdim)={1};
am = double(apc.MYDOUBLE);
am = am(rep{:});

else
 
Hm = mean(getPC(apc));
n=ndims(apc.MYDOUBLE);
nopcdim = setdiff(1:n,apc.pcdim);
s = size(apc.MYDOUBLE);
am = permute(double(apc.MYDOUBLE),[nopcdim , apc.pcdim]);
am = reshape(am,prod(s(1:end-1)),s(end))*Hm;
am = reshape(am,s(1:end-1));
  
end

