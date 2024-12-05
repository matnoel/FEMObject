function co = RFCORREL(varargin)
% function co = RFCORREL(varargin)
% See also EXPCORREL, EXP2CORREL
if isa(varargin{1},'RFCORREL')
    co=varargin{1};
else
    co.type = varargin{1} ;
    co.param = varargin{2} ;
    co = class(co,'RFCORREL');       
end
