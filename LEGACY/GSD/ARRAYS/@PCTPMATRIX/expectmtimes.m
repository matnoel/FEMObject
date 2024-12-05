function x = expectmtimes(x,varargin)

x = expectnodimmtimes([],x,varargin{:});
x = simplify(x);
