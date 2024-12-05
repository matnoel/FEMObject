function rep = iseval(ls)
% function rep = iseval(ls)

rep = 1;

for k=1:ls.n
    rep=rep & iseval(ls.LS{k});
end

