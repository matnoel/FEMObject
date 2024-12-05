function Hx = expand(H)

% Pense dans le cas scalaire uniquement (cf norme)
Hx=0; 
alpha = H.alpha;
for r=1:H.m
    exr=alpha(r);
    for d=1:H.dim
        exr = exr * expand(H.F{r,d});
    end
    Hx=Hx+exr;
end







