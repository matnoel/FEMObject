function x = icdf(u,P)

param=getparam(u);

B = (param.B);
A = (param.A);

P = min(1,max(0,P));

x = zeros(size(P));
rep = find(P==0);
x(rep) = A;

rep = find(P==1);
x(rep) = B;

rep = find(P<1 & P>0);
x(rep) = A*(B/A).^P(rep);

