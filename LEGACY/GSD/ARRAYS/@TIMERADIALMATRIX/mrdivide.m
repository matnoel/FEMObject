function w=mrdivide(u,v)

if isa(u,'PCRADIALMATRIX') & isa(v,'double')
 
    if all(size(v)==1)
    w=u;
    w.V = u.V/v;
    else
        
    error('mrdivide accepte uniquement un double de taille 1-by-1')    
    end
    
else
       error('mrdivide not defined') 
end

