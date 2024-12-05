function u = getmode(u,k)

if isa(u.value,'PCRADIALMATRIX')
    u.value=getmode(u.value,k);
else
    error('la PCTIMEMATRIX n''est pas RADIAL au niveau stochastique')
end