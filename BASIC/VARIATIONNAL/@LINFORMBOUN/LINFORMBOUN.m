function l = LINFORMBOUN(q,G,varargin)
% function l = LINFORMBOUN(q,G,k)
% l(v) = int_frontiere D^p(v) . k

% function l = LINFORM(q,G,k,pk)
% k est un champ defini aux noeuds
% l(v) = int_frontiere D^q(v) . D^pk(k)
% par defaut : pk=[], k est alors un simple parametre scalaire ou vectoriel


if nargin==0
    l=struct();
    l.G = [];
    l = class(l,'LINFORMBOUN',MULTILINFORM());
else
    if length(q)>1
        error('rentrer l''ordre');
    end
    l=struct();
    if nargin==1
        G=[];
    end
    l.G = G;
    m = MULTILINFORM(q,varargin{:});
    l = class(l,'LINFORMBOUN',m);
end

