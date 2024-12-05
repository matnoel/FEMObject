function s=num2str(u)

n=getnumber(RANDVARS(u));
s=[num2str(size(u,1)) '-by-' num2str(size(u,2)) ' ' class(u) ': Stochastic Dimensions = ' num2str(n{:}) ];