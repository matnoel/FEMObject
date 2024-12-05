function u = mean(u,varargin)
% function u = mean(u,varargin)

if ~israndom(u)
    for k=1:length(u.value)
        u.value{k} = mean(u.value{k},varargin{:});
    end
else
    for k=1:length(u.value)
        u.value{k} = MYDOUBLEND(mean(full(u.value{k})));
    end
end
