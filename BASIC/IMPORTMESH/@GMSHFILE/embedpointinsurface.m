function u = embedpointinsurface(u,numberpoint,numbersurface)
% function u = embedpointinsurface(u,numberpoint,numbersurface)

u = embed(u,'Point',numberpoint,'Surface',numbersurface);
