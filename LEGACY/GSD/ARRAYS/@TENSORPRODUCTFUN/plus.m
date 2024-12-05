function c = plus(a,b)

if isa(a,'TENSORPRODUCTFUN') && isa(b,'TENSORPRODUCTFUN')
    V = union(getproductspace(a),getproductspace(b));
    c = TENSORPRODUCTFUNSUM(V,a,b);
elseif isa(a,'double')
    c = plus(TENSORPRODUCTFUN(PRODUCTSPACE(),a),b);
elseif isa(b,'double')
    c = plus(a,TENSORPRODUCTFUN(PRODUCTSPACE(),b));
else
    error('pas programme')
end

