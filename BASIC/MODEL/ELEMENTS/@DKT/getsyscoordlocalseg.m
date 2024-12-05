function sys = getsyscoordlocalseg(elem,p)
%function sys = getsyscoordlocalseg(elem,p)
% obtenir les parametres des elements p

param = getparam(elem) ;
sys = param.syscoordlocalseg;

if nargin==2
    sys = getsubmatrix(sys,p,3);
end
