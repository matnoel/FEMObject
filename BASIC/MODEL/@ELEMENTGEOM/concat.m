function elems=concat(elem1,elem2)
% on rajoute des elements a un groupe d'element, si la connectivite est la
% meme l'element n'est pas ajoute.
elems=elem1;
elems.connec=[elems.connec;elem2.connec];
elems.numelem=[elems.numelem;elem2.numelem];

elems.syscoordlocal=concat(elems.syscoordlocal,elem2.syscoordlocal);
elems.syscoord=concat(elems.syscoord,elem2.syscoord);
elems = concatparam(elems,elem2);
elems.nbelem = size(elems.connec,1);