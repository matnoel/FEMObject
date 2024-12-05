function x = std(x,varargin)

Ex = expect(x);
x = expectnodimtimes([],x,x,varargin{:});
x = sqrt(simplify(x)-Ex.^2);

