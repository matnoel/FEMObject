function l = TRILINFORM(p,q,r,varargin)
% function a = TRILINFORM(p,q,k)
% a(u,v,w) = int k D^p(u) D^q(v) D^r(w)
% k : parametre
%
% function a = TRILINFORM(p,q,r,k,pk)
% k est un champ defini aux noeuds
% a(u,v) = int D^pk(k)  D^p(u) D^q(v) D^r(w)
% par defaut : pk=[], k est alors un simple parametre scalaire ou vectoriel

if nargin==0
    l = struct();
    l = class(l,'TRILINFORM',MULTILINFORM());

else
%     if isempty(varargin)
%     varargin={1 0};
%     end
    l = struct();
    m = MULTILINFORM([p,q,r],varargin{:});
    l = class(l,'TRILINFORM',m);
end
