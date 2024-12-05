function u = truncate(u,varargin)

if isa(u.value,'PCRADIALMATRIX')
    u.value=truncate(u.value,varargin{:});
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end