function a=double(apc,choix)

if isa(apc.value,'double') || isa(apc.value,'logical')
a = double(apc.value) ;

if nargin==2 && strcmp(choix,'reshape')
a = reshape(full(a),[apc.s size(apc.value,2)]);   
end

elseif isa(apc.value,'cell')
if numel(apc.value)==1
    a = apc.value{1};
else
    error('stockage cell de la multimatrix : utiliser cell')
end

else
    error('value of multimatrix must be double or cell ')

end
