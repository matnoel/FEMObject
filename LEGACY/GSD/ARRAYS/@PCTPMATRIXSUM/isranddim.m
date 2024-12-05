function r = isranddim(x,varargin)

r = isranddim(x.funs{1},varargin{:});
    
for k=2:getm(x)
  r = r | isranddim(x.funs{k},varargin{:});
end

