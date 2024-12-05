function [Fxy,x,y]=ecdf2D(Xs,Ys,nx,ny);
% function [Fxy,x,y]=ecdf2D(Xs,Ys,nx,ny)
% Xs et Ys, echantillons des variables X et Y
% nx, ny : nombre de points selon les dimensions x et y
% ou points d'evaluation de la CDF
% 
% x,y : grille d'evaluation de la CDF
% Fxy, valeurs de F(x,y)

Q = numel(Xs);
dx =  (-min(Xs)+max(Xs))/Q;
dy =  (-min(Ys)+max(Ys))/Q;
if length(nx)==1
    x = linspace(min(Xs)-dx,max(Xs)+dx,nx);
else
    x=nx;
    nx=length(x);
end
if length(nx)==1
    y = linspace(min(Ys)-dy,max(Ys)+dy,ny);
else
    y=ny;
    ny=length(y);
end

[x,y] = meshgrid(x,y);
x=x';
y=y';

Xs = repmat(Xs(:),1,nx*ny);
Ys = repmat(Ys(:),1,nx*ny);
xt = repmat(x(:)',Q,1);
yt = repmat(y(:)',Q,1);

Fxy = sum(Xs<=xt & Ys<=yt,1)/Q;
Fxy = reshape(Fxy,nx,ny);




