function rep = isalsmaterial(ls)

if ~isemty(ls.material) && isa(ls.material,'MATERIAL') 
    rep=1;
else
    rep=0;
end