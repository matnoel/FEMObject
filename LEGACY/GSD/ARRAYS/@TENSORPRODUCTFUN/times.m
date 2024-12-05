function c = times(a,b)

if isa(a,'TENSORPRODUCTFUN') && isa(b,'TENSORPRODUCTFUN')
    V = union(a.PRODUCTSPACE,b.PRODUCTSPACE);    
    [sa,repa] = ismember(V,a.PRODUCTSPACE);
    [sb,repb] = ismember(V,b.PRODUCTSPACE);
    
    c=cell(1,getnbdim(V));
    for i=1:getnbdim(V)
        if sa(i) && sb(i)
           c{i} = times(a.phi{repa(i)},b.phi{repb(i)}); 
        elseif sa(i) && ~sb(i)
           c{i} = a.phi{repa(i)};
        elseif ~sa(i) && sb(i)
           c{i} = b.phi{repb(i)};
        end
    end
    c = TENSORPRODUCTFUN(V,a.factor*b.factor,c{:});
else
    
    error('pas programme')
end

