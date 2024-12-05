function [N,rep] = getpointnextto(N,P)
% function [N,rep] = getpointnextto(N,P)

if isa(P,'POINT')
    if numel(P)>1
        error('rentrer un seul point')
    end
    d = double(distance(N,P));
    [d,rep] = min(d(:));
    N = getpoint(N,rep);
    
else
    error('pas programme')
end
