function X = random(rv,varargin)
% function X = random(rv,varargin)

for i=1:length(rv.Y)
    rv.Y{i} = random(rv.Y{i},varargin{:});
end
for i=1:length(rv.funparam)
    if isa(rv.funparam{i},'function_handle') || isa(rv.funparam{i},'inline')
        rv.funparam{i} = rv.funparam{i}(rv.Y{:});
    end
end

X = rv.X(rv.funparam{:});
X = random(X,varargin{:});
