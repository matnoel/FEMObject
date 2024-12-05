function gauss = calc_gauss(a,elem,order)

if nargin==2
   order = getorder(a.q,elem) + getorder(a.pk,elem);
end
    gauss = calc_gauss(elem,order);
    
function o = getorder(p,elem)
if isempty(p)
    o=0;
else
switch p
    case 0
       o = orderN(elem);
    case 1
       o = orderDN(elem); 
end
end


return


