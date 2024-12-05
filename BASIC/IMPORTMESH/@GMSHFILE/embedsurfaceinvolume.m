function u = embedsurfaceinvolume(u,numbersurface,numbervolume)
% function u = embedsurfaceinvolume(u,numbersurface,numbervolume)

u = embed(u,'Surface',numbersurface,'Volume',numbervolume);
