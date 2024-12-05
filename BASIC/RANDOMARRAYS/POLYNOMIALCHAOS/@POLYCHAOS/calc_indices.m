function PC = calc_indices(PC)
% function PC=calc_indices(PC)

PC.indices = PCbase_indices(PC.M,max(PC.n)-1,PC.typebase);
if length(PC.n)>=PC.M
    for k=1:PC.M
        rep=find(PC.indices(:,k)>PC.n(k)-1);
        PC.indices(rep,:)=[];
    end
else
    error('ordres mal definis')
end

