function u=orthogonalizeL(u,varargin)

if isa(u.value,'PCRADIALMATRIX')
    u.value=orthogonalizeL(u.value,varargin{:});
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end