function rvu = unique(rv)
% function rvu = unique(rv)

number = getnumber(rv);
[number,ia] = unique([number{:}]);
rvu=rv;
rvu.RV = rv.RV(ia);
rvu.M = length(rvu.RV);
rvu = sort(rvu);

[ok,ib] = ismember(rv,rvu);
ok=1;
for i=1:length(ib)
    ok = ok & (rv.RV{i}==rvu.RV{ib(i)});
end

if ok==0
    error('deux variables de meme numero sont differentes')
end

