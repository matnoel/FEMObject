function l = LINFORM(q,varargin)
% function l = LINFORM(q,k)
% l(v) = int_domaine D^q(v) . k
%
% function l = LINFORM(q,k,pk)
% k est un champ defini aux noeuds
% l(v) = int_domaine D^q(v) . D^pk(k)
% par defaut : pk=[], k est alors un simple parametre scalaire ou vectoriel

if nargin==0
    l=struct();
    l = class(l,'LINFORM',MULTILINFORM());
else
    if length(q)>1
        error('rentrer l''ordre');
    end
    l=struct();
    m = MULTILINFORM(q,varargin{:});
    l = class(l,'LINFORM',m);
end

