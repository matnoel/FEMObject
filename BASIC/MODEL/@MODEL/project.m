function S = project(S,varargin)
% function S = project(S,pc)
% on projette les LEVELSETS et les MATERIALS sur le POLYCHAOS pc
%
% See also MATERIALS/project, LEVELSETS/project, PCMATRIX/project, RANDVAR/project,
% RANDVARS/project

if israndom(S.ls);
    S.ls = project(S.ls,varargin{:});
end

mat = MATERIALS(S);
if israndom(mat)
    mat = project(mat,varargin{:});
    S = actualisematerials(S,mat);
end
