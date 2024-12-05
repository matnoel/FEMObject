function polys = unique(polys)
% function polys = unique(polys)

number = getnumber(polys);
[number,ia] = unique([number{:}]);
polys.h = polys.h(sort(ia));
polys.M = length(polys.h);

