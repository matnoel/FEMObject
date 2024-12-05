function elems = getelem(eleme,num,option,option2)
% function elems = getelem(eleme,num)
% num est la liste des numeros des elements (numerotation locale au groupe d'element)
% function elems=getelem(eleme,num,'global')
% num est la liste des numeros des elements (numerotation globale)


elems=eleme;

% [rep,ia,ib]=intersect(get(eleme,'numelem'),num);
if nargin==4 && strcmp('unique',option2)
    ia=unique(num);
else
    ia=num;
end



if nargin>=3 && strcmp('global',option)
    [temp,ia] = ismember(num,getnumber(eleme));
    ia = ia(find(temp));
end

ia=nonzeros(ia);


elems.connec = elems.connec(ia,:);
elems.numelem = elems.numelem(ia);

elems.syscoordlocal=getsyscoord(elems.syscoordlocal,ia);
elems.syscoord=getsyscoord(elems.syscoord,ia);

% elems.NODE = getnode(elems.NODE,unique(elems.connec));

elems = setparam(elems,getparam(elems,ia));
elems.nbelem = length(ia);


if ~isempty(getnumddl(elems))
    numddl = getnumddl(elems,ia);
    elems = setnumddl(elems,numddl);
end

