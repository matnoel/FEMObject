function ls=restrict(ls,varargin)

for i=1:ls.n
    ls.LS{i} = restrict(ls.LS{i},varargin{:});
end

