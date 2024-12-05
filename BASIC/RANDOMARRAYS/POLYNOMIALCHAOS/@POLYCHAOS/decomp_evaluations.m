function upc = decomp_evaluations(apc,gp,y)
% function u = decomp_evaluations(apc,gp,y)
% Decomposition of a function of random variables on polynomial chaos basis
% for given evaluations
% PC : POLYCHOAS or PCMODEL or FESTOMODEL or PCMATRIX or PCARRAY
% gp : gauss points
% y : matrix of evaluations of the random vector at gauss points
% y(i,j) : evaluations of the ith random variable at jth gauss point

PC = getPC(apc);

D = polyval(PC,gp.coord);
W = spdiags(gp.w(:),0,gp.nbgauss,gp.nbgauss);


upc = PCMATRIX(y*W*D,[size(y,1) 1],PC);

