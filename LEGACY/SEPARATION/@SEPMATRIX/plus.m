function w = plus(u,v)

w = u;
w.m = u.m+v.m;
w.F = [u.F ; v.F];
w.alpha = [u.alpha,v.alpha];

