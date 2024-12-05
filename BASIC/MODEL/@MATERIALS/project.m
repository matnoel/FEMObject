function mats = project(mats,varargin)
% function mats = project(mats,pc)
% projection de mats sur le chaos pc
% on applique project(.,pc) a tous les MATERIAL
% 
% See also MATERIAL/project

for k=1:mats.n
    mats.MAT{k} = project(mats.MAT{k},varargin{:});
end
