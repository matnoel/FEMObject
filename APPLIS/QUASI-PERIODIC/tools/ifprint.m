function [] = ifprint(condition, varargin)
% function [] = ifprint(condition, varargin)

if condition
    fprintf(varargin{:})
end

end