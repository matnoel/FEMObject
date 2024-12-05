function varargout=contourplot(ls,varargin)
%function varargout=contourplot(ls,varargin)

for i=1:length(ls.LS)
    contourplot(ls.LS{i},varargin{:});
end
