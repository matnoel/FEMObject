function w=getmatrix(u,k)

if nargin==1 && ( (isa(u.value,'double') && size(u.value,2)==1) || (isa(u.value,'cell') && numel(u.value)==1))
    if isa(u.value,'cell')
    w=reshape(u.value{1},u.s);
    else
    w=reshape(u.value,u.s);
    end
else
    
if length(k)>1
    error('getmatrix marche pour une seule matrice')
end
    

    if isa(u.value,'cell')
    w=reshape(u.value{k},u.s);
    else
    w=reshape(u.value(:,k),u.s);
    end
end