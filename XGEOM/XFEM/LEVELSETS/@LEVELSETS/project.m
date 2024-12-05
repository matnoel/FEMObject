function ls = project(ls,varargin)
% function ls = project(ls,pc)
% projection de ls sur le chaos pc
% on applique project(.,pc) a toutes les LEVELSET
% 
% See also LEVELSET/project

for k=1:ls.n
    ls.LS{k} = project(ls.LS{k},varargin{:});
end
