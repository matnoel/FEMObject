function elem=setstomesh(elem,h,e)
% function elem=setstomesh(elem,h)
% elem : groupe d'éléments
% h : tableau de cellules , h{e} : POLYFEND correspondant au maillage
% stochastique de l'élément e du groupe d'élément 'elem'
if nargin==2
elem.stomesh=h;
elseif nargin ==3
    elem.stomesh{e}=h;
else
    error('mauvais arguments')
end