function H = normalizefuns(H)
%  function H = normalizefuns(H)
%  Transfere la norme de toutes les HSM/SM 
%  dans les alpha correspondant.

alpha = H.alpha;
for r=1:H.m
    exr=1;
    for d=1:H.dim
        nrd      = norm(H.F{r,d});
        H.F{r,d} = (1/nrd) * H.F{r,d};
        exr      = exr * nrd;
        H.F{r,d} = normalizefuns(H.F{r,d});
    end
    alpha(r) = alpha(r)*exr;
end

H.alpha = alpha;