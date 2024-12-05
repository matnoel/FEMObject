function PC = projectdim(PC,dim,varargin)
% function PC = projectdim(PC,dim)
% Reduce PC to dimension dim

[ok,rep]=ismember(dim,PC);
if ~all(ok)
    error('Wrong dimensions')
end

PC.RANDPOLYS = PC.RANDPOLYS(rep);
d = setdiff(1:PC.M,dim);
ind = find(sum(PC.indices(:,d),2)>0);
PC.indices(ind,:) = [];
PC.indices = unique(PC.indices(:,rep),'rows','stable');
switch PC.typebase
    case 1
        PC.indices(:,end+1) = sum(PC.indices,2);
    case 2
        PC.indices(:,end+1) = max(PC.indices,[],2);
end
if ischarin('sort',varargin)
    PC.indices = sortrows(PC.indices,size(PC.indices,2):-1:1);
end
PC.M = size(PC.indices,2)-1;
PC.p = PC.p(rep);
PC.n = PC.n(rep);
PC.P = size(PC.indices,1)-1;
PC.masseuni = {};
PC.masse = {};
PC.metricuni = {};
PC.metric = {};