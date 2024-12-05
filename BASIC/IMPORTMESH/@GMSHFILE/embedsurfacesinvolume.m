function u = embedsurfacesinvolume(u,numbersurfaces,numbervolume)
% function u = embedsurfacesinvolume(u,numbersurfaces,numbervolume)

for k=1:length(numbersurfaces)
    u = embedsurfaceinvolume(u,numbersurfaces(k),numbervolume);
end
