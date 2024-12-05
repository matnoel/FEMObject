function M = changeelemnumber(M)
% function M = changeelemnumber(M)

M.nbelem=0;
groupgarde=zeros(1,M.nbgroupelem);
M.repelemingroupelem=zeros(0,2);
for j=1:M.nbgroupelem
    nbelem=getnbelem(M.groupelem{j});
    number = M.nbelem+[1:nbelem];
    M.groupelem{j}=setnumelem(M.groupelem{j},number);
    M.repelemingroupelem=[M.repelemingroupelem;...
        repmat(j,length(number),1),(1:length(number))'];
    
    if nbelem>0
        groupgarde(j)=1;
    end
    M.nbelem = M.nbelem + nbelem;
end

M.groupelem=M.groupelem(find(groupgarde));
M.nbgroupelem=length(M.groupelem);


