function [rep,pos] = isclassinvarargin(s,var)
% function [rep,pos] = isclassinvarargin(s,var)

warning('obsolete : replace by isclassin')
pos=[];
rep=0;
for i=1:length(var)
    if isa(var{i},s)
        rep=1;
        pos=[pos,i];
    end
end
