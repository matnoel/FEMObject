function l = BILINFORM(p,q,varargin)
% function a = BILINFORM(p,q,k)
% a(u,v) = int k D^p(u) D^q(v)
% k : parametre
%
% function a = BILINFORM(p,q,k,pk)
% k est un champ defini aux noeuds
% a(u,v) = int D^pk(k)  D^p(u) D^q(v)
% par defaut : pk=[], k est alors un simple parametre scalaire ou vectoriel

if nargin==0
    l=struct();
    l = class(l,'BILINFORM',MULTILINFORM());
else
    if length(q)>1
        error('rentrer l''ordre');
    end
    l=struct();
    m = MULTILINFORM([p,q],varargin{:});
    l = class(l,'BILINFORM',m);
end

