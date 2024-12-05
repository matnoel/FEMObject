function t = gettm(T)

t=T.t;
t = 1/2*(t(2:end)+t(1:end-1));

