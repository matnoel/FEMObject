function c = plus(a,b)

if isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUN')
    c=a;
    c.tensorfuns = [c.tensorfuns, {b}];  
    c.PRODUCTSPACE = union(a.PRODUCTSPACE,getproductspace(b));
    
elseif isa(a,'TENSORPRODUCTFUN') && isa(b,'TENSORPRODUCTFUNSUM')
    c=b;
    c.tensorfuns = [{a} , c.tensorfuns]; 

    c.PRODUCTSPACE = union(getproductspace(a),b.PRODUCTSPACE);
    
elseif isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUNSUM')
    c = a ; 
    c.tensorfuns = [a.tensorfuns , b.tensorfuns];
    c.PRODUCTSPACE = union(a.PRODUCTSPACE,b.PRODUCTSPACE);
elseif isa(a,'double') && numel(a)==1
    c = plus(TENSORPRODUCTFUN(0,a),b);
elseif isa(b,'double') && numel(b)==1
    c = plus(a,TENSORPRODUCTFUN(0,b));
else
    error('pas programme')
end

