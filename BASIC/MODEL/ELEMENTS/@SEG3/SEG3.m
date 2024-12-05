function elem = SEG3(node,numelem,connec,varargin)
% function elem = SEG3(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements
% option : 'BORD' si c'est un element de bord

if nargin==0    
    elemp = ELEMENTGEOM(1);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'SEG3',elemp,eleme);
else
    connec = reshape(connec,length(numelem),3);
    
    P1 = POINT(getnode(node,connec(:,1)));
    P2 = POINT(getnode(node,connec(:,2)));
    eX = P2 -P1 ;
    
    eX = normalize(eX);
    %param.L = distance(P1,P2);
    
    syscoordlocal = CARTESIAN1D(eX);
    
    elemp = ELEMENTGEOM(1,node,numelem,connec,syscoordlocal);
    eleme = ELEMENT(varargin{:});
    
    elem = struct();
    elem = class(elem,'SEG3',elemp,eleme);
end
