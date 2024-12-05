function mat = project(mat,varargin)
% function mat = project(mat,pc)
% projection de mat sur le chaos pc
% on applique project(.,pc) a tous les parametres
% 
% See also MATERIALS/project
 
mat.param = funrandomparam(mat.param,@project,varargin{:});
