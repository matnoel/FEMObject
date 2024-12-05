function DN = DNtri3(xnode)
% DN = DNtri3(xnode)

J = Jtri3(xnode);
DN = inv(J)*[-1,1,0;-1,0,1];
