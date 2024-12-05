function Xpc = eicdfpcproject2D(Xs,p,h1,h2,varargin)
% function Xpc = eicdfpcproject2D(Xs,p,h1,h2,'nbblock',N)
% Xs : echantillons de 2 VA : taille 2-by-N
% h1,h2 : RANDPOLY
% p : ordre du chaos
% nbblock : nombre de sous-ensemble pour l'estimation de la CDF conditinnelle
%      -> approx. de l'evenement {Y=y} par {Y in A}, avec A un sous-ensemble de l'ensemble des
%       realisations de Y

nbblock = getcharin('nbblock',varargin);
if isempty(nbblock)
    error('rentrer le nombre de blocks')
end

if nargin<=2 || isempty(h1) 
    h1 = POLYHERMITE();
end
if nargin<=2 || isempty(h2) 
    h2 = POLYHERMITE();
end
h = RANDPOLYS(h1,h2);
pc1 = POLYCHAOS(h1,p(1));
PC = POLYCHAOS(h,p);

x1pc = eicdfpcproject(Xs(1,:),p(1),h1);
x1pc = project(x1pc,PC);

Xs1 = Xs(1,:);
Xs2 = Xs(2,:);

ng = getcharin('nbgauss',varargin);
x2pc = decompfun(PC,ng,[],@(xi) myfun(xi,RANDVARS(PC),Xs1,Xs2,nbblock));

Xpc = full([x1pc;x2pc]);



function y = myfun(xi,RV,Xs1,Xs2,n)

d = (max(Xs1)-min(Xs1))/n;
X1 = eicdf(Xs1,cdf(RV{1},xi(:,1)));
y=zeros(size(xi,1),1);

for i=1:size(xi,1)
    rep = find(Xs1<=X1(i)+d/2 & Xs1>=X1(i)-d/2);
    y(i) = eicdf(Xs2(rep),cdf(RV{2},xi(i,2)));
end

return





