function x = calc_ximasse(x,varargin)

for i=1:length(x.funs)
   x.funs{i} = calc_ximasse(x.funs{i},varargin{:});
end

