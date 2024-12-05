function DN = DNtet4(xnode)
% DN = DNtet4(xnode)

J = Jtet4(xnode);
DN = inv(J)*[-1,1,0,0;-1,0,1,0;-1,0,0,1];
