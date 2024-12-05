function Z=getpccompo(apcr,p)

L = apcr.D * getpccompo(apcr.L,p);
if iscell(apcr.V)
Z = apcr.V{1}*L(1);
    for k=2:apcr.m
    Z = Z + apcr.V{i}*L(i); 
    end    
else
Z = double(apcr.V)*L;
Z = reshape(Z,size(apcr));
end