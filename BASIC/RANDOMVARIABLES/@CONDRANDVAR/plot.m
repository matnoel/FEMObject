function varargout = plot(u,varargin)
% function [P,x] = plot(u,courbe)
% trace la pdf ou cdf ou icdf de u (CONDRANDVAR)
% courbe  = 'pdf' ou 'icdf' ou 'cdf' ('pdf' par defaut)
% plot(u,courbe,'npts',n) trace la courbe avec n points (100 par defaut)
% plot(u,courbe,'bar','npts',n) trace des barres.
% plot(u,courbe,'nbs',N) nombre de tirage (10000 par defaut)
% tous les autres arguments sont passes en argument de plot ou bar
%
% See also RANDVAR/pdfplot, RANDVAR/cdfplot, RANDVAR/icdf, PCMATRIX/plot, sampleplot

if ~ischarin('multipdf',varargin)
    N=getcharin('nbs',varargin,1e5);
    varargin = delcharin('nbs',varargin);
    
    us = random(u,N,1);
    varargout=cell(1,nargout);
    [varargout{:}]=sampleplot(us,varargin{:})   ;
else
    if length(u.Y)>1
        error('affichage d''une pdf en 3D ou plus non programme')
    end
    
    ny = getcharin('npts',varargin,100);
    nx = ny;
    dy=getdomainborne(u.Y{1});
    y=linspace(dy(1),dy(2),ny);
    
    pY = pdf(u.Y{1},y);
    
    for i=1:length(u.funparam)
        if isa(u.funparam{i},'function_handle') | isa(u.funparam{i},'inline')
            u.funparam{i}=u.funparam{i}(y);
        end
    end
    X = u.X(u.funparam{:});
    dx=getdomainborne(X);
    x=linspace(dx(1),dx(2),nx);
    
    pXY = zeros(nx,ny);
    for i=1:nx
        pXY(i,:) = pdf(X,x(i));
    end
    pY = repmat(pY(:)',nx,1);
    pXY = pXY.*pY;
    
    [x,y] = meshgrid(y,x);
    
    
    varargin=delonlycharin({'pdf','cdf','icdf','bar','multipdf'},varargin);
    varargin=delcharin({'nbs','nbbar','npts'},varargin);
    
    surf(x,y,pXY,varargin{:})
    xlabel('x')
    ylabel('y')
    
    if nargout>=1
        varargout{1}=pXY;
    end
    if nargout>=2
        varargout{2}=xy;
    end
    
end