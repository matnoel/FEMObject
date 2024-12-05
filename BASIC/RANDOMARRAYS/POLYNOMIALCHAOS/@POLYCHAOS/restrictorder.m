function PC = restrictorder(PC,order,dim,typebase)
% function PC = restrictorder(PC,degree,dim,typebase)

if nargin<=2 || isempty(dim)
    dim = 1:PC.M;
end

[ok,rep]=ismember(dim,PC);
if ~all(ok)
    error('Wrong dimensions')
end

if nargin<=3
    typebase = 1;
end

for i=1:length(rep)
    if isa(getpoly(PC,rep(i)),'POLYFE')
        error('It does not work with POLYFE')
    end
end

switch typebase
    case 1
        if numel(order)>1
            error('For typebase = 1, enter a single polynomial degree')
        end
        ind = find(sum(PC.indices(:,rep),2)>order);
        PC.indices(ind,:) = [];
        PC.p(rep) = min(PC.p(rep),order);
        PC.n(rep) = min(PC.n(rep),order+1);
    case 2
        if numel(order)==1
            order = repmat(order,1,length(rep));
        end
        if numel(order)~=numel(rep)
            error('For typebase 2, enter as many polynomial degres as dimensions')
        end
        ind = [];
        for i=1:length(rep)
            ind = union(ind,find(PC.indices(:,rep(i))>order(i)));
        end
        PC.indices(ind,:) = [];
        PC.p(rep) = min(PC.p(rep),order);
        PC.n(rep) = min(PC.n(rep),order+1);
end

PC.P = size(PC.indices,1)-1;
PC.masseuni = {};
PC.masse = {};
PC.metricuni = {};
PC.metric = {};
