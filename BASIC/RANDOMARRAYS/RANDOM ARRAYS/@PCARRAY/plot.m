function [P,x] = plot(u,varargin)
% function [P,x] = plot(u,varargin)
% trace la pdf ou cdf ou icdf de u (PCARRAY)
% plot(u,courbe,varargin) avec courbe  = 'pdf' ou 'icdf' ou 'cdf' (pdf par defaut)
% pdfplot(u,arg1,...) arguments de plot
% pdfplot(u,'npts',n,arg1,...) trace la courbe avec n points (100 par defaut)
% pdfplot(u,'bar','npts',n,arg1,...) trace des barres. arg1 ... arguments de bar
% pdfplot(u,'nbs',N) nombre de tirage (10000 par defaut)

pos=1:nargin-1;

n=getcharin('npts',varargin,100);
N=getcharin('nbs',varargin,1e5);
tracebar=ischarin('bar',varargin);
nbbar=getcharin('nbbar',varargin,1);
varargin=delonlycharin('bar',varargin);
varargin = delcharin({'npts','nbs','nbbar'},varargin);

[rep,post]=ischarin({'pdf','cdf','icdf'},varargin);
courbe = varargin{post(find(rep))};
varargin = delonlycharin({'pdf','cdf','icdf'},varargin);

Nmax=1e4;
switch courbe
    case 'cdf'
        us=randomlimit(u,N,Nmax);
        [P,x] = ecdf(us);
    case 'pdf'
        us=randomlimit(u,N,Nmax);
        
        Nb = floor(1+3.3*log(N)/log(10))*nbbar ;
        %Nb = floor(sqrt(N)*nbbar);
        
        x=linspace(min(us),max(us),Nb+1);
        P = prob(x,us);
        x=[x(1:end-1)+x(2:end)]/2;
        
        %[P,x] = ksdensity(us);
    otherwise
        error(' ')
end

if tracebar
    bar(x,P,varargin{:})
else
    plot(x,P,varargin{:})
end

function [P] = prob(x,us)
P=zeros(1,length(x)-1);
for k=1:length(x)-1
    P(k)=length(find(us>= x(k) & us<x(k+1)))/(x(k+1)-x(k))/length(us) ;
end


function us = randomlimit(u,N,Nmax)
us=random(u,N-floor(N/Nmax)*Nmax);
for k=1:floor(N/Nmax)
    us=[us(:);reshape(random(u,Nmax),Nmax,1)];
end
