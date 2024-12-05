function x = expecttimes(x,varargin)
% function x = expecttimes(x,varargin)

x = expectnodimtimes([],x,varargin{:});
x = simplify(x);
