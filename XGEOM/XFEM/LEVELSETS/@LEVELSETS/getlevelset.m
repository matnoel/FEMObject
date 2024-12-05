function level = getlevelset(levels,k)
% function level = getlevelset(levels,k)

if nargin==1
    level = levels.LS;
else
    [ok,rep] = ismember(k,levels);
    if length(k)~=length(rep)  || any(rep==0)
        error('wrong number')
    end
    level = levels.LS(rep);
end

if length(level)==1
    level = level{1};
elseif isempty(level)
    level = [];
end