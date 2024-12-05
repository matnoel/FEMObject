function a=double(apc,choix)

a = apc.value ;

if nargin==2 & strcmp(choix,'reshape')
a = reshape(full(a),[apc.s size(apc.value,2)]);   
end
