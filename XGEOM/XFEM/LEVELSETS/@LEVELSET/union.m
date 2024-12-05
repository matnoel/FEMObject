function ls1 = union(ls1,ls2)
% function ls1 = union(ls1,ls2)

if iseval(ls1) && iseval(ls2)
    ls1.value = min(ls1.value,ls2.value);
else
    ls1 = LEVELSET(@union,ls1,ls2);
end
