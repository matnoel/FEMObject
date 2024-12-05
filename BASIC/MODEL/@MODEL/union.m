function M = union(M,M2,varargin)
% function M = union(M,M2,option)
% union de deux modeles M et M2 en eliminant les doubles noeuds (fonction unique)
% On ajoute a M les elements de M2 (fonction addelem)
%
% si option = 'duplicate' : pas d'elimination des doubles noeuds
% si option = 'norenum' : pas de renumerotation des noeuds
% sinon, on considere que les modeles ont leur numerotation propre
% et union cree alors un modele avec renumerotation des noeuds
%
% See also MODEL/intersect, MODEL/setdiff, MODEL/unique, MODEL/addelem

if ~isa(M,'MODEL') || ~isa(M2,'MODEL')
    error('rentrer deux MODEL pour effectuer union')
end

options = delclassin('MODEL',varargin);
M = addelem(M,M2,options{:});

[rep,pos] = isclassin('MODEL',varargin);
if rep
    M = union(M,varargin{pos(1)},varargin{setdiff(1:length(varargin),pos(1))});
end
