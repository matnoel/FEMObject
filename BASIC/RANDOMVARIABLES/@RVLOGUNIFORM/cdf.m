function P = cdf(u,x)

param=getparam(u);

B = (param.B);
A = (param.A);

P = zeros(size(x));
rep = find(x>=B);
P(rep)=1;
rep = find(x>=A & x<B);
P(rep)=(log(x(rep))-log(A))/(log(B)-log(A));

