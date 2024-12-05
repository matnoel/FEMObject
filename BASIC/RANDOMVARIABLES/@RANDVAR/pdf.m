function P = pdf(u,x)

P = pdf(u.type,x,u.param{:,2});
