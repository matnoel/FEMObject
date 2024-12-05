function ls=restrict(ls,M1,M2)

if iseval(ls)
posnode = getpos(getnode(M1),getnumber(getnode(M2)));
ls.value = ls.value(posnode,:);
end

