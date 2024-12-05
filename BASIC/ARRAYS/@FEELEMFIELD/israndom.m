function ok = israndom(u)
% function ok = israndom(u)

ok = 1;
for k=1:length(u.value)
    ok = ok & israndom(u.value{k});
end
