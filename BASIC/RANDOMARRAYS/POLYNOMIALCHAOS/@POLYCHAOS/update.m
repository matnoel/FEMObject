function PC = update(PC)
% function PC = update(PC)
% Update POLYCHAOS PC after any modification of its multi-index set

PC.M = size(PC.indices,2)-1;
PC.p = max(PC.indices(:,1:PC.M),[],1);
% PC.n = PC.p+1;
for k=1:PC.M
    PC.n(k) = length(unique(PC.indices(:,k)));
end
PC.P = size(PC.indices,1)-1;
PC.typebase = PC.typebase;
PC.masseuni = {};
PC.masse = {};
PC.metricuni = {};
PC.metric = {};
