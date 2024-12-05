function a = PCMATRIX(apc)

s = size(apc.MYDOUBLE);
s(apc.pcdim)=[];
if length(s)==1
    s=[s,1];
end
value=double(apc);
rep = 1:ndims(apc.MYDOUBLE);
rep(apc.pcdim) = ndims(apc.MYDOUBLE);
rep(end) = apc.pcdim;
value = permute(value,rep);
value = reshape(value,prod(s),length(apc.POLYCHAOS));
a = PCMATRIX(value,s,apc.POLYCHAOS);
