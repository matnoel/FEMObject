function c = mylocnode(c,elem,node,N)
c=full(c);
if isempty(c)
    c=1;
elseif numel(c)~=getdim(elem)
    if size(c,1)==getnbnode(node)
        c=N*localize(elem,c);  
    else
        c = mylocelem(c,elem);  
        if size(c,1)==getdim(elem)^2
            c = reshape(c,[getdim(elem),getdim(elem),sizeND(c)]);    
        end
    end
    if size(c,1)==1 
        c = c';
    end
end
