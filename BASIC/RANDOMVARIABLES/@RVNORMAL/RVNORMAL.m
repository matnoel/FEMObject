function rv = RVNORMAL(varargin)
% function rv = RVNORMAL(mu,sigma)

rv=struct();
switch nargin
case 0
    param.mu = 0;
    param.sigma = 1;
otherwise
    param.mu = varargin{1};
    param.sigma = varargin{2};
end


rvp = RANDVAR('normal',param);
rv = class(rv,'RVNORMAL',rvp);

