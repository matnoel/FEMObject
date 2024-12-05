function u=normalizeL(u,varargin)

if isa(u.value,'PCRADIALMATRIX')
    u.value=normalizeL(u.value,varargin{:});
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end