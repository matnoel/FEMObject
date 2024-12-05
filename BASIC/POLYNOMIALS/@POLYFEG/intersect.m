function H = intersect(H1,H2)
% function H = intersect(H1,H2)

if ~(isa(H1,'POLYFEG') && isa(H2,'POLYFEG'))
    error('pas programme')
elseif ~(H1.rv == H2.rv)
    error('germes differents')
else
    H.POLYFE = intersect(H1.POLYFE,H2.POLYFE);
end

