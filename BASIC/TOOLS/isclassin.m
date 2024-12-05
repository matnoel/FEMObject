function [rep,pos] = isclassin(s,var)
% function [rep,pos] = isclassin(s,var)

pos=[];
rep=0;
for i=1:length(var)
    if isa(var{i},s)
        rep=1;
        pos=[pos,i];
    end
end
