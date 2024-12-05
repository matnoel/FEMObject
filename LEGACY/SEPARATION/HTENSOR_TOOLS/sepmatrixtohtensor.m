function [H] = sepmatrixtohtensor(S,varargin)
% function [H] = sepmatrixtohtensor(S,varargin)
%   Automatically vectorize an operator
%   varargin: options passed to the htensor constructor

if size(S.F{1},2)>1
    isoper = 1;
else
    isoper = 0;
end


if isoper
    for i=1:S.m
        for j=1:S.dim
            S.F{i,j}=S.F{i,j}(:);
        end
    end
end

H=gathervectors(normalizefactors(S));
H=htensor(H.F',varargin{:});

end
