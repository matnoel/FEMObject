function varargout = quiver(u,M,ampl,varargin)
% function varargout = quiver(u,M,ampl,varargin)
% affichage d'un champ de vecteur defini aux noeuds
% u : FENODEFIELD
% M : MODEL
% ampl : amplification du champ de vecteur

coord = getcoord(M.node);
vec = double(u.value);
dimvec = size(vec,2);
dim = M.dim;
coord = [coord,zeros(size(coord,1),dimvec-dim)];

if nargin==2
    ampl = {};
else
    ampl = {ampl};
end
rep = find(sum(abs(vec),2)~=0);

switch dimvec
    case 2
      H = quiver(coord(rep,1),coord(rep,2),vec(rep,1),vec(rep,2),ampl{:},varargin{:});  
    case 3
      H = quiver3(coord(rep,1),coord(rep,2),coord(rep,3),vec(rep,1),vec(rep,2),vec(rep,3),ampl{:},varargin{:});    
end

if nargout>=1
    varargout{1} = H;
end
  