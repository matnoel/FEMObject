function mat = ELAS_BEAM(varargin)
% function mat = ELAS_BEAM('parameter',value)
% Parameters:
% 'E' Young modulus
% 'NU' Poisson ratio
% 'IY' Planar second moment of area along Y axis
% 'IZ' Planar second moment of area along Z axis
% 'IX' Polar second moment of area along X axis
% 'RHO' mass density
% 'S' cross-section area
% 'MU' damping factor

mat = struct();

param = struct(varargin{:});

matp = MATERIAL('ELAS_BEAM',param);
mat = class(mat,'ELAS_BEAM',matp);
