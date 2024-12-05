function u = random(u,varargin)
% function u = random(u,varargin)

if ~israndom(u)
    error('FEELEMFIELD is not random')
end

for k=1:length(u.value)
    u.value{k} = random(full(u.value{k}),varargin{:});
    if ~isa(u.value{k},'MULTIMATRIX')
        u.value{k} = MYDOUBLEND(u.value{k});
    end
end
