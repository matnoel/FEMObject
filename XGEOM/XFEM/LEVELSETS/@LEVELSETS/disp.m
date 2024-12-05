function disp(levels,varargin)
% function disp(levels,varargin)

for k=1:levels.n
    if ~isempty(levels.LS{k})
        fprintf('{%d} ',k)
        disp(levels.LS{k},varargin{:})
    end
end
