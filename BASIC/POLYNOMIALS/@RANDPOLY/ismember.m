function [rep,i] = ismember(k,u)
% function [rep,i] = ismember(k,u)

if ~isa(k,'double')
    k = getnumber(k);
    if isa(k,'cell')
        k = [k{:}];
    end
end
if ~isa(u,'double')
    u = getnumber(u);
    if isa(u,'cell')
        u = [u{:}];
    end
end

[rep,i] = ismember(k,u);


