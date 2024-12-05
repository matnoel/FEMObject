function [ok,rep] = israndom(ls)
% function [ok,rep] = israndom(ls)

rep = zeros(1,ls.n);

for k=1:ls.n
    rep(k) = israndom(ls.LS{k});
end
ok = any(rep);
rep = find(rep);

