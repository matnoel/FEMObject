function mat = ELAS_SHELL(varargin)
% function mat = ELAS_SHELL('parameter',value)
% Parameters:
% 'E' Young modulus
% 'NU' Poisson ratio
% 'RHO' mass density
% 'DIM3' thickness
% 'k' shear correction factor
% 'MU' damping factor

mat = struct();

param = struct(varargin{:});

if ~isfield(param,'k')
    param.k = 5/6; % shear correction factor
end

matp = MATERIAL('ELAS_SHELL',param);
mat = class(mat,'ELAS_SHELL',matp);
