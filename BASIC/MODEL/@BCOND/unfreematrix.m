function v = unfreematrix(BC,u)
% function v = unfreematrix(BC,u)

v = unfreevector(BC,u,1);
v = unfreevector(BC,v,2);
