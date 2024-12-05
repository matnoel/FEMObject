function elem = CUB8(node,numelem,connec,varargin)
% function elem = CUB8(node,numelem,connec,'material',mat)
% node : objet de type NODE contenant les noeuds de l'element
% numelem : numero des elements

if nargin==0
    elemp = ELEMENTGEOM(3);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'CUB8',elemp,eleme);
elseif nargin==1 && isa(node,'NODE') && numel(node)==8
    elem = CUB8(node,1,1:8);  
elseif nargin==1 && isa(node,'double') && size(node,1)==8
    elem = CUB8(NODE(node),1,1:8);  
else
    syscoordlocal = CARTESIAN3D();
    syscoordlocal = repmat(syscoordlocal,[1,1,size(connec,1)]);
    syscoord = syscoordlocal;
    
    elemp = ELEMENTGEOM(3,node,numelem,connec,syscoordlocal,syscoord);
    eleme = ELEMENT(varargin{:});
    
    elem = struct();
    elem = class(elem,'CUB8',elemp,eleme);
end
