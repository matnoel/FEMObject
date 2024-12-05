function n=nbvect(H)

n=0;
for r=1:H.m
    nr=0;
    for d=1:H.dim
        nrd=nbvect(H.F{r,d});
        nr=nr+nrd;
    end
    n=n+nr;
end