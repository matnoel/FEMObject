function v = variance(u)
% function v = variance(u)

if isconditional(u)
    error('moyenne pas calulee pour variables conditionnelles')
else
    [m,v]=rvstat(u);
end