function ls1=setdiff(ls1,ls2)
if iseval(ls1) && iseval(ls2)
    ls1.value = max(ls1.value,-ls2.value);
else
    ls1 = LEVELSET(@setdiff,ls1,ls2);
end
