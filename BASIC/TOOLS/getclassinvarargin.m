function [v,num] = getclassinvarargin(s,var,default)
% function [v,num] = getclassinvarargin(s,var,default)

warning('obsolete : replace by getclassin')

pos=[];
rep=0;
for i=1:length(var)
    if isa(var{i},s)
        rep=1;
        pos=[pos,i];
    end
end
if length(pos)==1
    v=var{pos};
elseif length(pos)>1
    v=var(pos);
else
    v=[];
end

num = length(pos);

if isempty(pos) && nargin>2
    v=default;
end
