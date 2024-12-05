function a = PCMATRIX(a)

if isradial(a)
    warning('passage en PC du stockage du PCMYDOUBLEND')
    a = expand(a);
end

a = PCMATRIX(MULTIMATRIX(a.V,a.stodim),a.L);