function m = getm(u)

if isa(u.value,'PCRADIALMATRIX')
    m=getm(u.value);
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end