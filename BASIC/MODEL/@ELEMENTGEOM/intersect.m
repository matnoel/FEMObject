function [elem,node] = intersect(elem1,elem2,node)
% function [elem,node] = intersect(elem1,elem2,node)

if isa(elem1,class(elem2))
    loc =  ismember(sort(elem1.connec,2),sort(elem2.connec,2),'rows')    ;
    elem=getelem(elem1,find(loc));
    if nargin==3 && nargout==2
        node = getnode(node,getnumnode(elem));
    end
else
    error('les elements ne sont pas du meme type')
end