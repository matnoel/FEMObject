function x = normalize(x)

if ~isa(x,'double') || length(size(x))>2
    error('pas prevu')
else
    x = x/norm(x);   
end
