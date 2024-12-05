function a = randomeval(x,xi,RV)
% function a = randomeval(x,xi,RV)

if nargin==2
    Hs = randomeval(x.POLYCHAOSTP,xi);
else
    Hs = randomeval(x.POLYCHAOSTP,xi,RV);
end

a = randomevalwithpoly(x,Hs);


