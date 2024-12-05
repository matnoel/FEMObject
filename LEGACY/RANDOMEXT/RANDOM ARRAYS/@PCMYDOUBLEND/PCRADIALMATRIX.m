function a = PCRADIALMATRIX(a)

if ~isradial(a)
    warning('passage en RADIAL du stockage du PCMYDOUBLEND')
    a = convertradial(a);
end

a = PCRADIALMATRIX(MULTIMATRIX(a.V,a.stodim),[],a.L);