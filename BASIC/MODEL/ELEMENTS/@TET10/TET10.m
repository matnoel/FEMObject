function elem = TET10(node,numelem,connec,varargin)
% function elem = TET10(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements
% option : 'BORD' si c'est un element de bord
%          'DEFO' 'CONT'

if nargin==0
    elemp = ELEMENTGEOM(3);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'TET10',elemp,eleme);
elseif nargin==1 && isa(node,'NODE') && numel(node)==10
    elem = TET10(node,1,1:10);  
elseif nargin==1 && isa(node,'double') && size(node,1)==10
    elem = TET10(NODE(node),1,1:10);  
else
    syscoordlocal = CARTESIAN3D();
    syscoord = syscoordlocal;
    
    elemp = ELEMENTGEOM(3,node,numelem,connec,syscoordlocal,syscoord);
    eleme = ELEMENT(varargin{:});
    
    elem = struct();
    elem = class(elem,'TET10',elemp,eleme);
end
