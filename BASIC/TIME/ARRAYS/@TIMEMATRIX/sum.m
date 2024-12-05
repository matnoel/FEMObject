function v = sum(u,varargin)

if isa(u.value,'cell')
    v=u;
    for i=1:length(v.value)
        v.value{i}=sum(v.value{i},varargin{:});
    end
    v.s = size(v.value{1});
else
v=u;
U = MULTIMATRIX(double(u.value),u.s,[length(u.TIMEMODEL),prod(sizem(u.value))]);
v.value = sum(U,varargin{:});
v.s = size(v.value);
v.value = double(v.value);
if isa(u.value,'MULTIMATRIX')
v.value = MULTIMATRIX(v.value,[prod(v.s),length(u.TIMEMODEL)],sizem(u.value));   
end
end
    