function u=sum(u,k)
if any(k==u.multidim)
    error('pas de sommation sur la dimension multi')
end
u.value=sum(u.value,k);