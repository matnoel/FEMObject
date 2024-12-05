function rep = isparam(mat,i)
% function rep = isparam(mat,i)

if isa(i,'char')
    [rep,pos] = ischarin(i,mat.param(:,1));
else
    error('rentrer un char')
end


