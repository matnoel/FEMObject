function levels = setlevelset(levels,ls,k)
% function levels = setlevelset(levels,ls,k)

if nargin==2
    levels.LS = LEVELSETS(ls);
else
    ls= setnumber(ls,k);
    levels.LS{k} = ls;
end

levels.n = length(levels.LS);
