function gauss = trianglequad(order)
% function gauss = trianglequad(order)

n = ceil((order+1)/2);
[x,w] = simplexquad(n,[0,0;1,0;0,1]);

gauss.coord = x;
gauss.w = w;
gauss.nbgauss = length(w);
