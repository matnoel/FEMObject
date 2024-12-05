function varargout = vectorplot(u,varargin)
% function varargout = vectorplot(u,varargin)

for l=1:nargin-1
    switch class(varargin{l})
        case 'POINT'
            v = varargin{l};
        case 'char'
            plotoptions = varargin{l};
        case 'double'
            ampl = varargin{l};
    end
    if ~exist('plotoptions','var')
        plotoptions = 'b';
    end
    if ~exist('v','var')
        v = POINT([0,0]);
    end
    if ~exist('ampl','var')
        ampl = 1;
    end
end

u = squeeze(double(u));
v = squeeze(double(v));

H = quiver(v(1,:),v(2,:),u(1,:),u(2,:),ampl,plotoptions);

if nargout>=1
    varargout{1} = H;
end
