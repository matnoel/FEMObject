function s = calc_opmat(S,varargin)
% function s = calc_opmat(S,varargin)
% S : MODEL

if ischarin('node',varargin)
    s = calc_elemfield(S,@opmatnode,varargin{:});
else
    s = calc_elemfield(S,@opmat,varargin{:});
end
