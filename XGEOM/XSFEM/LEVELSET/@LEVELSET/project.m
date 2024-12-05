function ls = project(ls,varargin)
% function ls = project(ls,pc)
% on projette la LEVELSET sur le POLYCHAOS pc
% - si ls est une levelset et ls.value est une PCMATRIX
%   on appelle project(ls.value,pc)
% - sinon on applique project(.,pc) a tous les parametres de la
%   levelset
%
% See also LEVELSETS/project, PCMATRIX/project, RANDVAR/project,
% RANDVARS/project


if isalevelset(ls)
    ls.value = project(ls.value,varargin{:});
else
    if isa(ls.value{1,1},'function_handle')
        levels = ls.value{1,2};
        for k=1:length(levels)
            levels{k} = project(levels{k},varargin{:});
        end
        ls.value{1,2} = levels ;
    else
        ls.value = funrandomparam(ls.value,@project,varargin{:});
    end
end
