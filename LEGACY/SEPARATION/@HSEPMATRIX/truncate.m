function u = truncate(u,i)
u.alpha = u.alpha(i);
u.F = u.F(i,:);
u.m = length(u.alpha);


