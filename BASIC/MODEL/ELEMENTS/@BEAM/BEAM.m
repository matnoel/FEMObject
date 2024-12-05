function elem = BEAM(node,numelem,connec,varargin)
% function elem = BEAM(node,numelem,connec,'material',mat,'param',eY)
%
% node : objet de type NODE contenant les noeuds de l'element
% matnum : numero du materiau
% connec : table de connectivite
% numelem : numero des elements
% eY (uniquement en 3D) : eX etant le vecteur tangent a l'element, 
%       eY est un vecteur permettant de definir le plan principal (eX,eY) de la poutre

if nargin==0
    elemp = SEG2();
else
    elemp = SEG2(node,numelem,connec,varargin{:});
end
elem = struct();
elem = class(elem,'BEAM',elemp);
elem = setlocal(elem);
