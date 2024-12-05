function m = mean(u)
% function m = mean(u)

if isconditional(u)
    error('moyenne pas calulee pour variables conditionnelles')
else
    m=rvstat(u);
end