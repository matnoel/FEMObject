function mats = unique(mats)
% function mats = unique(mats)

number = getnumber(mats);
[number,ia] = unique([number{:}]);
mats.MAT = mats.MAT(sort(ia));
mats.n = length(mats.MAT);

