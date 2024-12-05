function elem = QUA4(node,numelem,connec,varargin)
% function elem = QUA4(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% numelem : numero des elements
% connec : table de connectivite
% option : 'BORD' si c'est un element de bord
%          'DEFO' 'CONT'
if nargin==0
    elemp = ELEMENTGEOM(2);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'QUA4',elemp,eleme);
elseif nargin==1 && isa(node,'NODE') && numel(node)==4
    elem = QUA4(node,1,1:4);
elseif nargin==1 && isa(node,'double') && size(node,1)==4
    elem = QUA4(NODE(node),1,1:4);
else
    P1 = POINT(getnode(node,connec(:,1)));
    P2 = POINT(getnode(node,connec(:,2)));
    P4 = POINT(getnode(node,connec(:,4)));
    eX = normalize(P2-P1);
    eY = normalize(P4-P1);
    if getindim(node)==3
        eZ = cross(eX,eY);
        eZ = normalize(eZ);
        eY = cross(eZ,eX);
        syscoordlocal = CARTESIAN2D(eX,eY);
        syscoord = CARTESIAN3D(eX,eY,eZ);
    else
        eY = rot2D(eX,pi/2);
        syscoordlocal = CARTESIAN2D(eX,eY);
        syscoord = syscoordlocal;
    end
    
    elemp = ELEMENTGEOM(2,node,numelem,connec,syscoordlocal,syscoord);
    eleme = ELEMENT(varargin{:});
    
    elem = struct();
    elem = class(elem,'QUA4',elemp,eleme);
end
