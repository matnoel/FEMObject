function u = normalize(u,varargin)
u.alpha = u.alpha/norm(u,varargin{:});
