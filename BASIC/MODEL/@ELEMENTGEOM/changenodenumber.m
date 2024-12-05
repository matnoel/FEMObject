function elem = changenodenumber(elem,oldnumber,newnumber)
% function elem = changenodenumber(elem,oldnumber,newnumber)

[a,newconnec] = ismember(elem.connec,oldnumber);
if ~all(a)
    error('la connectivite est associee a des noeuds inexistants')
end
elem.connec = reshape(newnumber(newconnec(a)),size(elem.connec));
elem.nbelem = size(elem.connec,1);
