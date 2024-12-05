function ls = actualise(ls,varargin)
% function ls = actualise(ls,varargin)

for i=1:length(ls.LS)
    ls.LS{i} = actualise(ls.LS{i},varargin{:});
end
