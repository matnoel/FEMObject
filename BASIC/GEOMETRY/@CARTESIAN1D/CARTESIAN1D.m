function s = CARTESIAN1D(varargin)
% function s = CARTESIAN1D(varargin)

if nargin==0
    s.name = 'Coordonnees cartesiennes 1D';
    sp = SYSCOORD(1);
    s = class(s,'CARTESIAN1D',sp);
elseif nargin==1 && isa(varargin{1},'CARTESIAN1D')
    s = varargin{1};
else
    s.name = 'Coordonnees cartesiennes 1D';
    base = getcompo(VECTEUR(varargin{1}));
    axis = {'X'};
    sp = SYSCOORD(base,axis);
    s = class(s,'CARTESIAN1D',sp);
end
