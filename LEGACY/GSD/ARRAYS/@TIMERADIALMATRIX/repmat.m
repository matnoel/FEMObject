function u = repmat(u,varargin)

switch nargin
    case {3,2}
u.V = repmat(u.V,varargin{:});

    otherwise
        
        error('pas le bon nombres d''arguments')
end


