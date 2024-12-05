function varargout = plot(u,varargin)
% function varargout = plot(u,varargin)

numnode=ischarin('numnode',varargin);
varargin = delonlycharin('numnode',varargin);
if numnode
    color = getcharin('color',varargin,'b');
    plotnumber(u,'color',color);  
end
color = delcharin('color',varargin);

H = plot(u.POINT,varargin{:});

if nargout>=1
    varargout{1}=H;
end


