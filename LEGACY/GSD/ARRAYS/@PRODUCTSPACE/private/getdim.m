function a = getdim(a)

if isa(a,'PRODUCTSPACE')
    a=a.dim;
end
