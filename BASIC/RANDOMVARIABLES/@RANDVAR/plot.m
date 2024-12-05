function varargout = plot(u,varargin)
% function [P,x] = plot(u,courbe)
% trace la pdf ou cdf ou icdf de u (RANDVAR)
% courbe  = 'pdf' ou 'icdf' ou 'cdf' ('pdf' par defaut)
% plot(u,courbe,'npts',n) trace la courbe avec n points (100 par defaut)
% plot(u,courbe,'bar','npts',n) trace des barres.
% plot(u,courbe,'nbs',N) nombre de tirage (10000 par defaut)
% tous les autres arguments sont passes en argument de plot ou bar
%
% See also RANDVAR/pdfplot, RANDVAR/cdfplot, RANDVAR/icdf, PCMATRIX/plot, sampleplot

d=getdomainborne(u);
n = getcharin('npts',varargin,100);

if ischarin('cdf',varargin)
    x=linspace(d(1),d(2),n);
    P=cdf(u,x);
elseif ischarin('icdf',varargin)
    x=linspace(0,1,n);
    P=icdf(u,x);
else
    x=linspace(d(1),d(2),n);
    P = pdf(u,x);
end

tracebar=ischarin('bar',varargin);
varargin=delonlycharin({'pdf','cdf','icdf','bar'},varargin);
varargin=delcharin({'npts'},varargin);

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

