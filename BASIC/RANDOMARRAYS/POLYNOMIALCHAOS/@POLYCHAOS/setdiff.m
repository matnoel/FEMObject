function [PC,i1] = setdiff(PC1,PC2,varargin)
% function PC = setdiff(PC1,PC2)
% soustraction d'un POLYCHAOS PC2 a un POLYCHAOS PC1
%
% See also POLYCHAOS/intersect, POLYCHAOS/union

PC.typebase = getcharin('typebase',varargin,min(PC1.typebase,PC2.typebase));

if ~polycmp(PC1.RANDPOLYS,PC2.RANDPOLYS)
    error('Different polynomial bases')
end

PC = PC1;
ind1 = getindices(PC1);
ind2 = getindices(PC2);
[PC.indices,i1]=setdiff(ind1(:,1:PC.M),ind2(:,1:PC.M),'rows','stable');
switch PC.typebase
    case 1
        PC.indices(:,PC.M+1) = sum(PC.indices(:,1:PC.M),2);
    case 2
        PC.indices(:,PC.M+1) = max(PC.indices(:,1:PC.M),[],2);
end
if ischarin('sort',varargin)
    PC.indices = sortrows(PC.indices,size(PC.indices,2):-1:1);
end
PC.p = max(PC.indices(:,1:PC.M),[],1);
for k=1:PC.M
    PC.n(k) = length(unique(PC.indices(:,k)));
end
PC.P = size(PC.indices,1)-1;
PC.masseuni = {};
PC.masse = {};
PC.metricuni = {};
PC.metric = {};
