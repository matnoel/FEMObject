function ls1 = intersect(ls1,ls2)
% function ls1 = intersect(ls1,ls2)

if iseval(ls1) && iseval(ls2)
    ls1.value = max(ls1.value,ls2.value);
else
    ls1 = LEVELSET(@intersect,ls1,ls2);
end
