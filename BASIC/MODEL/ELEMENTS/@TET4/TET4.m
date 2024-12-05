function elem = TET4(node,numelem,connec,varargin)
% function elem = TET4(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements
% option : 'BORD' si c'est un element de bord
%          'DEFO' 'CONT'

if nargin==0
    elemp = ELEMENTGEOM(3);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'TET4',elemp,eleme);
elseif nargin==1 && isa(node,'NODE') && numel(node)==4
    elem = TET4(node,1,1:4);  
elseif nargin==1 && isa(node,'double') && size(node,1)==4
    elem = TET4(NODE(node),1,1:4);  
else
    syscoordlocal = CARTESIAN3D();
    syscoord = syscoordlocal;
    
    elemp = ELEMENTGEOM(3,node,numelem,connec,syscoordlocal,syscoord);
    eleme = ELEMENT(varargin{:});
    
    elem = struct();
    elem = class(elem,'TET4',elemp,eleme);
end
