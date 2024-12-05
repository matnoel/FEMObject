function u = max(u,varargin)


if isa(u.value,'cell')
for i=1:length(u.value)
u.value{i} = max(u.value{i},varargin{:}); 
end
elseif isa(u.value,'FEELEMFIELD')  

    if nargin==3 & varargin{2}==2
    warning('ambiguite sur le max')
    end

u.value = max(u.value,varargin{:});


if isa(u.value,'FEELEMFIELD')
warning('ambiguite')
u.s = [size(u.value,1),1];
else
u.s = [size(u.value,1),1];
u.value = reshape(u.value,prod(u.s),length(u.TIMEMODEL));  
end

elseif isa(u.value,'double')

if all(u.s==1)
if nargin==2
    u.value = max(u.value,varargin{:});
else
    u = max(u.value,varargin{:});
end
else
v = MULTIMATRIX(u.value,u.s,[length(u.TIMEMODEL),1]); 
v = max(v,varargin{:});
u.s = size(v);
u.value = double(v);
end

else
error('pas programme')
end
    

