function elem=getgroupelem(M,p)

if nargin==1
    elem = M.groupelem;
else
    elem = M.groupelem{p};
end

