function [ok,rep] = israndom(mats)
% function [ok,rep] = israndom(mats)

rep = zeros(1,mats.n);
for k=1:mats.n
    rep(k) = israndom(mats.MAT{k});
end
ok = any(rep);
rep = find(rep);

