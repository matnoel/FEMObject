function PC = restrictdim(PC,dim,varargin)
% function PC = restrictdim(PC,dim)
% Restrict PC to dimension dim

[ok,rep]=ismember(dim,PC);
if ~all(ok)
    error('Wrong dimensions')
end

PC.RANDPOLYS = PC.RANDPOLYS(rep);
PC.indices = unique(PC.indices(:,rep),'rows','first');
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
end
% H = PC.RANDPOLYS(rep) ;
% p = PC.p(rep);
% 
% PC = POLYCHAOS(H,p,'typebase',PC.typebase);