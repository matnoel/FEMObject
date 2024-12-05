function [gauss,order] = calc_gauss(a,elem,order)

if nargin==2
    order=0;
    for i=1:length(a.p)
        order = order + getorder(a.p(i),elem);
    end
    order = order + getorder(a.pk,elem);
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


