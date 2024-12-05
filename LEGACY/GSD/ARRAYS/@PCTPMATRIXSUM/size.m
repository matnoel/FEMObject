function s = size(x,varargin)

if getm(x)>0
s = size(x.funs{1},varargin{:});
else
s = [];    
end