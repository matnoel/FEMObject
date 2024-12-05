function u=uniqueL(u,varargin)

if isa(u.value,'PCRADIALMATRIX')
    u.value=uniqueL(u.value,varargin{:});
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end