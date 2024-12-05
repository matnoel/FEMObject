function l = MIXEDBILINFORM(p,q,dir,varargin)
% function a = MIXEDBILINFORM(p,q,k)
% a(u,v) = int k D^p(u) D^q(v)
% k : parametre
% u, v definis sur deux MODEL differents
%
% function a = BILINFORM(p,q,k,pk)
% k est un champ defini aux noeuds
% a(u,v) = int D^pk(k)  D^p(u) D^q(v)
% par defaut : pk=[], k est alors un simple parametre scalaire ou vectoriel

if nargin==0
    l=struct();
    l.dir = [];
    l = class(l,'MIXEDBILINFORM',MIXEDMULTILINFORM());
else
    if length(q)>1
        error('rentrer l''ordre');
    end
    l=struct();
    if nargin>=3
        l.dir = dir;
    else
        l.dir=[];
    end
    m = MIXEDMULTILINFORM([p,q],varargin{:});
    l = class(l,'MIXEDBILINFORM',m);
end

