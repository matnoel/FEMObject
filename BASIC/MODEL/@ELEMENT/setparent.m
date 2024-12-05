function elem = setparent(elem,parent)
% function elem = setparent(elem,parent)
% 
if isa(parent,'ELEMENT')
elem.param.parent = class(parent);
elseif isa(parent,'char')
elem.param.parent = parent;    
else
error(' ')
end





