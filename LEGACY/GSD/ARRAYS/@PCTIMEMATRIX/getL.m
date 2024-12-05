function L = getL(u,varargin)

if isa(u.value,'PCRADIALMATRIX')
    L=getL(u.value,varargin{:});
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end