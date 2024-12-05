function varargout = fullsubplotbordure(n,m,k,bord,factx,facty,varargin)
% function fullsubplot(n,m,k,bord,factx,facty,varargin)
% fonction permettant d'obtenir des subplot d'une taille souhait�e
% n,m,k : 3 arguments classiques de subplot
% factx et facty des reels compris entre 0 et 1 indiquant la taille
% de chaque sous-figure (en x et y)

if isempty(m)
    ntemp = ceil(sqrt(n));
    m = ceil(n/ntemp);
    n = ntemp;
end

if nargin<5 || isempty(factx)
    factx=0.8;
end
if nargin<6 || isempty(facty)
    facty=factx;
end

j=1+mod(-1+k,m);
i=(k-j)/m+1;

margx=(1-factx)/2/m;
margy=(1-facty)/2/n;

checkfact(factx);
checkfact(facty);
checkfact(margx);
checkfact(margy);

sizex = factx/(m);
sizey = facty/n;

posx = (1-bord)/2 + bord*(margx*(2*j-1)+(j-1)*sizex);
posy = 1-(1-bord)/2-bord*(margy*(2*i-1)+(i)*sizey);

H = subplot('position',[posx,posy,sizex,sizey],varargin{:});
if nargout==1
    varargout{1} = H ;  
end

function checkfact(fact)

if fact<0 || fact>1
    error('le facteur doit etre entre 0 et 1')
end

return
