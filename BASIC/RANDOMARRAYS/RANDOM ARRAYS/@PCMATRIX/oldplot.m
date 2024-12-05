function varargout = plot(u,varargin)
% function [P,x] = plot(u,courbe)
% trace la pdf ou cdf ou icdf de u (PCMATRIX) 
% courbe  = 'pdf' ou 'icdf' ou 'cdf' ou 'rs' ('rs' par defaut)
% 'rs' trace la surface de reponse en fonction des variables de base du chaos
% plot(u,courbe,'npts',n) trace la courbe avec n points (100 par defaut)
% plot(u,courbe,'bar','npts',n) trace des barres. 
% plot(u,courbe,'nbs',N) nombre de tirage (10000 par defaut)    
% tous les autres arguments sont passes en argument de plot ou bar
%
% See also sampleplot, PCMATRIX/pdfplot, PCMATRIX/cdfplot, PCMATRIX/rsplot, RANDVAR/plot



pos=1:nargin-1;

n=getcharin('npts',varargin,100);
N=getcharin('nbs',varargin,1e5);
tracebar=ischarin('bar',varargin);
nbbar=getcharin('nbbar',varargin,1);
varargin=delonlycharin('bar',varargin);
varargin = delcharin({'npts','nbs','nbbar'},varargin);

[rep,post]=ischarin({'pdf','cdf','icdf','rs'},varargin);
if ~any(rep)
 courbe = 'rs'    ;
else
    courbe = varargin{post(find(rep))};
end

varargin = delonlycharin({'pdf','cdf','icdf','rs'},varargin);

switch courbe
    case 'cdf'
us=randomlimit(u,N);
[P,x] = ecdf(us);
    case 'pdf'
us=randomlimit(u,N);
Nb = floor(1+3.3*log(N)/log(10))*3*nbbar ;
%Nb = floor(sqrt(N)*nbbar);

x=linspace(min(us)-10*eps,max(us)+10*eps,Nb+1);
P = prob(x,us);
x=[x(1:end-1)+x(2:end)]/2;
%[P,x] = ksdensity(us);
    case 'icdf'
        error('pas programme')
    case 'rs'
RV = RANDVARS(u);
if length(RV)==1
    x = getdomainborne(RV{1});
    x=linspace(x(1),x(2),n);
    P = randomeval(u,x(:));
else
    error('surface de reponse pas programmee en 2D')
end

end

if tracebar
    bar(x,P,varargin{:})    
else
plot(x,P,varargin{:})   
end

if nargout>=1
    varargout{1}=P;
end
if nargout>=2
    varargout{2}=x;
end


function [P] = prob(x,us,Nmax)
P=zeros(1,length(x)-1);
for k=1:length(x)-1
P(k)=length(find(us>= x(k) & us<x(k+1)))/(x(k+1)-x(k))/length(us) ;  
end
return

function us = randomlimit(u,N,Nmax)
us=zeros(N,1);
for k=1:floor(N/Nmax)
us((k-1)*Nmax+[1:Nmax])=random(u,Nmax,1);
end
us(floor(N/Nmax)*Nmax+1:end)=random(u,N-floor(N/Nmax)*Nmax,1);
return

