function Z = getmatrixatstep(a,i)

if isa(a.value,'cell')
    Z=a.value{i};
elseif israndom(a.value)
    Z=reshape(a.value(:,i),a.s);
else
    error('pas prevu')
end