function mat = ELAS_SHELL_ISOT_TRANS(varargin)
% function mat = ELAS_SHELL_ISOT_TRANS('parameter',value)
% Parameters:
% 'ET' transverse Young modulus
% 'NUT' transverse Poisson ratio
% 'GL' longitudinal shear modulus
% 'RHO' mass density
% 'DIM3' thickness
% 'k' shear correction factor
% 'MU' damping factor
% The axis of rotational symmetry (axis L) corresponds to the axis normal to the plate mid-surface (axis Z)
% The plane of isotropy (any axis T) corresponds to the plate mid-surface (axes X and Y)

mat = struct();

param = struct(varargin{:});

if ~isfield(param,'k')
    param.k = 5/6; % shear correction factor
end

matp = MATERIAL('ELAS_SHELL_ISOT_TRANS',param);
mat = class(mat,'ELAS_SHELL_ISOT_TRANS',matp);
