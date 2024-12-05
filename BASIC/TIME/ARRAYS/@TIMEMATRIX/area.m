function area(u,varargin)

[t,rep]=gettplot(u);

switch class(u.value)
    case 'double'
area(t,u.value(:,rep),varargin{:});  
    case 'cell'
if ~isa(u.value{1},'double') && ~size(u.value{1},2)==1
    error('pas programme')
end
u.value = [u.value{:}];
area(u,varargin{:});
    otherwise
        class(u.value)
        error('pas programme')

end
