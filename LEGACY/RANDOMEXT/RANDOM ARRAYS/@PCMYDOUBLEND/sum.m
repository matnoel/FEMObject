function u=sum(u,k)
if any(k==u.stodim)
    error('pas de sommation sur la dimension stochastique')
end
u.V=sum(u.V,k);