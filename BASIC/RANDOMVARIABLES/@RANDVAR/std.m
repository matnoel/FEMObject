function s = std(u)
if isconditional(u)
error('moyenne pas calulee pour variables conditionnelles')
else
[m,v]=rvstat(u);
s=sqrt(v);
end