function M=convertelem(M,elemtype,varargin)
%function M=convertelem(M,elemtype,varargin)
% convertit les elements en elements elemtype
% M : MODEL
% elemtype : nom de l'element que l'on veut 
% varargin : options pour la defininition de l'element (voir les options de addelem)
%
% function M=convertelem(M,elemtype,'selgroup',selgroup)
% selgroup : numero des groupes d'elements a convertir

selgroup = getcharin('selgroup',varargin,1:M.nbgroupelem);

for j=selgroup
M.groupelem{j} = convert(M.groupelem{j},M.node,elemtype,varargin{:});
end

M = changeelemnumber(M);
