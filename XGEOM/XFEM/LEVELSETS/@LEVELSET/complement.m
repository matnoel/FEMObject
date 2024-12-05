function ls=complement(ls)
if iseval(ls)
    ls.value = -ls.value ;
else
    ls = LEVELSET(@complement,ls);
end
