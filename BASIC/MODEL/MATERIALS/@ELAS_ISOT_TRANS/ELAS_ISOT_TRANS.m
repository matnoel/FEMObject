function mat = ELAS_ISOT_TRANS(varargin)
% function mat = ELAS_ISOT_TRANS('parametre',valeur)
% Parameters:
% 'AXISL' axis of rotational symmetry (axis L)
% 'AXIST' any axis in the plane of isotropy (axis T)
% 'EL' longitudinal Young modulus
% 'ET' transverse Young modulus
% 'NUL' longitudinal Poisson ratio
% 'NUT' transverse Poisson ratio
% 'GL' longitudinal shear modulus
% 'RHO' mass density
% 'DIM3' thickness in the 2D case
% 'S' cross-section area in the 1D case
% 'MU' damping factor

mat = struct();

param = struct(varargin{:});

if ~isfield(param,'DIM3')
    param.DIM3 = 1;
end

matp = MATERIAL('ELAS_ISOT_TRANS',param);
mat = class(mat,'ELAS_ISOT_TRANS',matp);

