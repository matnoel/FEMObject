function c = prodscal(a,b,varargin)

if isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUNSUM')
    c = TENSORPRODUCTSUM();
    for i=1:length(a.tensorfuns)
    for j=1:length(b.tensorfuns)
    c = c + prodscal(a.tensorfuns{i},b.tensorfuns{j},varargin{:});    
    end
    end
    c = simplify(c);
   
    
elseif isa(a,'TENSORPRODUCTFUN') && isa(b,'TENSORPRODUCTFUNSUM')
    c = b;   
    for i=1:length(b.tensorfuns)
    c.tensorfuns{i} = prodscal(a,b.tensorfuns{i},varargin{:});    
    end
    c = simplify(c);
   
elseif isa(a,'TENSORPRODUCTFUNSUM') && isa(b,'TENSORPRODUCTFUN')   
    
    c = a;
    for i=1:length(a.tensorfuns)
    c.tensorfuns{i} = prodscal(a.tensorfuns{i},b,varargin{:});    
    end
    c = simplify(c);
    
else
    
    error('prodscal non defini')
    
end
