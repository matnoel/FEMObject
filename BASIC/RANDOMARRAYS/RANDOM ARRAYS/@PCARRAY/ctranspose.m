function apc=ctranspose(apc)

n=ndims(apc.MYDOUBLE);
switch n
    case 2
nopcdim = apc.pcdim;
apc.pcdim = mod(nopcdim,2)+1;
apc.MYDOUBLE = apc.MYDOUBLE';
    case 3
repnopc=setdiff(1:3,apc.pcdim);
rep=1:3;
rep(repnopc(1))=repnopc(2);
rep(repnopc(2))=repnopc(1);
apc.MYDOUBLE = permute(apc.MYDOUBLE,rep);        
    otherwise
        error('transpose not defined')
end
