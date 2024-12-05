function elems=addconnec(elem1,elem2)
% on rajoute des elements a un groupe d'element, si la connectivite est la
% meme l'element n'est pas ajoute.
elems=elem1;
elems.connec=[elems.connec;elem2.connec];
elems.connec=unique(elems.connec,'rows');
elems.nbelem = size(elems.connec,1);