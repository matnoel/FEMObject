function x = rdivide(a,b)

if isa(a,'PCTPMATRIX') && isa(b,'double')
    x = a;
    x.phi0 = rdivide(x.phi0,b);
    
else
    error('pas prevu')
    
end
