function hm=mean(h,liste)

param = get(h,'param');
hm = sparse(length(liste),1);
rep = find(liste<param.n);

I = mod(liste(rep),param.n)+1;
dx = param.I(I,2)-param.I(I,1);
hm(rep) = sqrt(dx);