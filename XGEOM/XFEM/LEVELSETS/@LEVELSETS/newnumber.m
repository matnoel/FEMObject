function n = newnumber(u,exclu)
% function n = newnumber(u,exclu)

if nargin==1
    exclu = {};
else
    if ~isa(exclu,'cell')
        exclu = num2cell(exclu{:});
    end
end
number = [getnumber(u),exclu];
allnumber = [number{:}];
n = min(setdiff(1:length(number)+1,allnumber));

