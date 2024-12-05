function a = randomeval(x,xi,RV)
% function a = randomeval(x,xi,RV)

if nargin==2
    Hs = randomeval(x.POLYCHAOSTP,xi);
else
    Hs = randomeval(x.POLYCHAOSTP,xi,RV);
end

a = randomevalwithpoly(x.funs{1},Hs);
for i=2:length(x.funs)
    a = a + randomevalwithpoly(x.funs{i},Hs);
end
