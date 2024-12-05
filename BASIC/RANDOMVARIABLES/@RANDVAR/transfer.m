function bx = transfer(a,b,ax)
% function bx = transfer(a,b,ax)
% a et b : RANDVAR    
% bx=F_b^-1(F_a(ax))
% ax,bx : n-by-1 double
% 
% See also RANDVARS/transfer

if ~isa(b,'RANDVAR')
    error('second argument doit etre une RANDVAR')
end

bx = icdf(b,cdf(a,ax));

