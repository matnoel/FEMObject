function c = times(a,b)

if isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUN')
    c=a;
    for i=1:a.M
    c.tensorfuns{i} = times(a.tensorfuns{i},b);  
    c.PRODUCTSPACE = union(c.PRODUCTSPACE,getproductspace(c.tensorfuns{i}));
    end
    
elseif isa(a,'TENSORPRODUCTFUN') && isa(b,'TENSORPRODUCTFUNSUM')
    c=b;
    for i=1:b.M
    c.tensorfuns{i} = times(a,b.tensorfuns{i});  
    c.PRODUCTSPACE = union(c.PRODUCTSPACE,getproductspace(c.tensorfuns{i}));
    end
elseif isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUNSUM')
    c = TENSORPRODUCTFUNSUM();
    c.tensorfuns = cell(1,a.M*b.M);
    for i=1:a.M
        for j=1:b.M
            k = (i-1)*b.M+j;
    c.tensorfuns{k} = times(a.tensorfuns{i},b.tensorfuns{j});  
    c.PRODUCTSPACE = union(c.PRODUCTSPACE,getproductspace(c.tensorfuns{k}));
        end
    end
  
else   
    error('pas programme')
end

