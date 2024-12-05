function x = expect(x,varargin)
% function x = expect(x,varargin)

x = expectnodim([],x,varargin{:});
x = simplify(x);
