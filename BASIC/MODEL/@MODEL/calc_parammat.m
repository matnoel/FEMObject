function s = calc_parammat(S,varargin)
% function s = calc_parammat(S,varargin)
% S : MODEL

if ischarin('node',varargin)
    s = calc_elemfield(S,@parammatnode,varargin{:});
else
    s = calc_elemfield(S,@parammat,varargin{:});
end
