function varargout = plot(u,varargin)

varargout=cell(1,nargout);
[t,rep]=gettplot(u);

switch class(u.value)
    case 'double'
        [varargout{:}]=plot(t,u.value(:,rep),varargin{:});
    case 'cell'
        if ~isa(u.value{1},'double') && ~size(u.value{1},2)==1
            error('pas programme')
        end
        u.value = [u.value{:}];
        [varargout{:}]=plot(u,varargin{:});
    otherwise
        error('pas programme')
        
end
