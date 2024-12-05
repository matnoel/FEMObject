function c = mtimes(a,b)

if isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUN')
    c=a;
    for i=1:length(a.tensorfuns)
    c.tensorfuns{i} = mtimes(a.tensorfuns{i},b);  
    c.PRODUCTSPACE = union(c.PRODUCTSPACE,getproductspace(c.tensorfuns{i}));
    end
    
elseif isa(a,'TENSORPRODUCTFUN') && isa(b,'TENSORPRODUCTFUNSUM')
    c=b;
    for i=1:length(b.tensorfuns)
    c.tensorfuns{i} = mtimes(a,b.tensorfuns{i});  
    c.PRODUCTSPACE = union(c.PRODUCTSPACE,getproductspace(c.tensorfuns{i}));
    end
elseif isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUNSUM')
    c = TENSORPRODUCTFUNSUM();
    c.tensorfuns = cell(1,length(a.tensorfuns)*length(b.tensorfuns));
    for i=1:length(a.tensorfuns)
        for j=1:length(b.tensorfuns)
            k = (i-1)*length(b.tensorfuns)+j;
    c.tensorfuns{k} = mtimes(a.tensorfuns{i},b.tensorfuns{j});
    c.PRODUCTSPACE = union(c.PRODUCTSPACE,getproductspace(c.tensorfuns{k}));
        end
    end
  
else   
    error('pas programme')
end

