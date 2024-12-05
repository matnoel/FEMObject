function u = getmydoublend(u,k)

    rep = cell(1,length(size(u.value)));
    rep(:)={':'};
    rep{u.multidim}=k;

if length(k)>1
    u = u.value(rep{:});
else
    u.value = u.value(rep{:});
end