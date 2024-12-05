function V = getV(u,varargin)

if isa(u.value,'PCRADIALMATRIX')
    V=getV(u.value,varargin{:});
    V = TIMEMATRIX(V,u.TIMEMODEL,size(u));
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end