function x = mrdivide(a,b)

if isa(a,'PCTPMATRIX') && isa(b,'double')
    x = a;
    x.phi0 = mrdivide(x.phi0,b);
    
else
    error('pas prevu')
    
end
