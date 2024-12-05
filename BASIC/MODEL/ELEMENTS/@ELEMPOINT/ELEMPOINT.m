function elem = ELEMPOINT(node,numelem,connec,varargin)
% function elem = ELEMPOINT(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite (vecteur de taille n-by-1)
% numelem : numero des elements
% option : 'BORD' si c'est un element de bord

if nargin==0
    elemp = ELEMENTGEOM(0);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'ELEMPOINT',elemp,eleme);
elseif nargin==1 && isa(node,'NODE')
    numelem = getnumber(node);
    elem = ELEMPOINT(node,getnbnode(node),numelem(:));  
else
    syscoordlocal = getsyscoord(node);
    syscoord = syscoordlocal;
    
    elemp = ELEMENTGEOM(0,node,numelem,connec,syscoordlocal,syscoord);
    eleme = ELEMENT(varargin{:});
    
    elem = struct();
    elem = class(elem,'ELEMPOINT',elemp,eleme);
end
