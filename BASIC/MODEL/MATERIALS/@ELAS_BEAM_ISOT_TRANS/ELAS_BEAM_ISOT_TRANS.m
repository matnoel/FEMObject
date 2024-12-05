function mat = ELAS_BEAM_ISOT_TRANS(varargin)
% function mat = ELAS_BEAM_ISOT_TRANS(varargin)
% Parameters:
% 'EL' longitudinal Young modulus
% 'NUT' transverse Poisson ratio
% 'GL' longitudinal shear modulus
% 'IY' Planar second moment of area along Y axis
% 'IZ' Planar second moment of area along Z axis
% 'IX' Polar second moment of area along X axis
% 'RHO' mass density
% 'S' cross-section area
% 'MU' damping factor
% The axis of rotational symmetry (axis L) corresponds to the beam axis (axis X)

mat = struct();

param = struct(varargin{:});

matp = MATERIAL('ELAS_BEAM_ISOT_TRANS',param);
mat = class(mat,'ELAS_BEAM_ISOT_TRANS',matp);
