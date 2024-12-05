function [se,B] = sigmaexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)
% function [se,B] = sigmaexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)

warning('faux')
N = calc_N(elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
u = N*qe;
gradu = B*qe;
k = evalparampc(mat,'k',[],elem,xnode,xgauss);

se = 2*k*gradu-(1-u).*u;

if ~israndom(se)
    se = se*expect(l1,l2);
else
    se = expect(se,l1,l2);
end
