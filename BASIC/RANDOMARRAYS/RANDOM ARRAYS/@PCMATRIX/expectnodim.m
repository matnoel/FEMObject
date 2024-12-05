function uk = expectnodim(k,u)
% function uk = expectnodim(k,u)

PC = getPC(u);
PCk = restrictdim(PC,k);
ind = getindices(PC);

sdim = setdiff(1:getM(PC),k);

if isempty(sdim)
    uk = u;
    return
end


for j=sdim
    H = mean(restrictdim(PC,j));
    H = H(ind(:,j)+1);
    if j==sdim(1)
        PSI = H;
    else
        PSI = PSI.*H;
    end
end

v = double(u);
uk = sparse(size(v,1),getP(PCk)+1);

if length(k)==1
    for i=0:getP(PCk)
        I = find(ind(:,k)==i);
        uk(:,i+1)=v(:,I)*PSI(I);
    end
    
else
    indk=getindices(PCk);
    for i=1:size(indk,1)
        I = ismember(ind(:,k),indk(i,1:end-1),'rows');
        uk(:,i)=v(:,I)*PSI(I);
    end
end

uk = PCMATRIX(uk,size(u),PCk);



