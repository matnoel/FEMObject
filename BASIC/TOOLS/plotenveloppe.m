function fillhandle = plotenveloppe(t,Y,face,varargin)
% function plotenveloppe(t,Y,facecolor)
% Y : matrice dont les colonnes sont des echantillons d'une fonction de t
% t : abscisse
% facecolor : couleur de remplissage
% 
% function plotenveloppe(Y,facecolor,'edgecolor',edgecolor)
% edgecolor : couleurs des contrours ('k' par defaut)
%
% function plotenveloppe(Y,facecolor,'facealpha',alpha,'edgealpha','alphaedge')
% alpha et edgealpha : coefficients de transparence (0.5 par defaut)

if nargin<=2
    face = 'y';
end

if strcmp(class(Y),'MULTIMATRIX')
    Y = double(Y);
elseif ~isa(Y,'double')
    error('pas programme')
end

if any(size(Y)==1)
    Y = reshape(Y,1,numel(Y));
    Ymax=max(Y,0);
    Ymin=-max(-Y,0);
else
    Ymin = min(Y,[],2)';
    Ymax = max(Y,[],2)';
end


edge = getcharin('edgecolor',varargin,'k');
transp = getcharin('facealpha',varargin,0.5);
transpedge = getcharin('edgealpha',varargin,transp);

options={};
if ischarin('linewidth',varargin)
    options=[options , {'linewidth',getcharin('linewidth',varargin)}];
end

t=t(:)';
t=[t,fliplr(t)];
Y = [Ymin,fliplr(Ymax)];
fillhandle = fill(t,Y,face);

if transp==1 && transpedge==1
    set(fillhandle,'EdgeColor',edge,options{:});
else
    set(fillhandle,'EdgeColor',edge,'FaceAlpha',transp,'EdgeAlpha',transpedge,options{:});    
end
