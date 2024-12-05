function Xpc = eicdfpcproject2D(Xs,p,h1,h2,varargin)
% function Xpc = eicdfpcproject2D(Xs,p,h1,h2,'nbblock',N)
% Xs : echantillons de 2 VA : taille 2-by-N
% h1,h2 : RANDPOLY
% p : ordre du chaos
% nbblock : nombre de sous-ensemble pour l'estimation de la CDF conditinnelle
%      -> approx. de l'evenement {Y=y} par {Y in A}, avec A un sous-ensemble de l'ensemble des
%       realisations de Y

N = size(Xs,2);
nbblock = getcharin('nbblock',varargin,floor((1+3.3*log(N)/log(10))*2));


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

%d = (max(Xs1)-min(Xs1))/n;

X1 = eicdf(Xs1,cdf(RV{1},xi(:,1)));

y=zeros(size(xi,1),1);

[F12,x1,x2] = ecdf2D(Xs1,Xs2,n,n);
x1m = x1(2:end,1)/2+x1(1:end-1,1)/2;


% figure(453)
% clf
% surf(x1,x2,F12);

D1F21 = (diff(F12,[],1)./(x1(2:end,:)-x1(1:end-1,:)));
%pX1 = repmat(ksdensity(Xs1,x1m),1,size(x1,2));
pX1 = repmat(D1F21(:,end),1,size(x1,2));
F2vert1 = D1F21./pX1;

% figure(454)
% clf
% surf(repmat(x1m,1,size(x1,2)),x2(2:end,:),F2vert1);

for i=1:length(y)
    rep = find(x1m(:,1)>=X1(i),1);
    if isempty(rep) || rep==1
        rep = 2;
    end

    F = F2vert1(rep,:);

    F = min(F,1);
    F = max(F,0);

%figure(23)
%clf
%plot(x2(1,:),F,'.');
%hold on
%plot(x2(1,:),F,'.');
%pause(1)

    y(i) = invF(x2(1,:),F,cdf(RV{2},xi(i,2))) ;
end
return


function x = invF(X,F,y)

if any(y>1 | y<0)
    error('la icdf n''est definie que sur [0,1]')
end
x=zeros(size(y));


for i=1:length(y)

    r1 = find(F>=y(i),1);    
    if isempty(r1)
        x(i) = X(end);
    else
        x(i) = X(r1);
    end

end

return



