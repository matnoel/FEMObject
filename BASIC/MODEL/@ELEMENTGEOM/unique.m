function elem = unique(elem)
% function elem = unique(elem)

[connecu,a,b] = unique(sort(getconnec(elem),2),'rows');
if size(connecu,1)~=getnbelem(elem)
    elem = getelem(elem,a);
end