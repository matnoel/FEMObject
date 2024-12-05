function levels = unique(levels)
% function levels = unique(levels)

number = getnumber(levels);
[number,ia] = unique([number{:}]);
levels.LS = levels.LS(sort(ia));
levels.n = length(levels.LS);

