function u = embedpointsinsurface(u,numberpoints,numbersurface)
% function u = embedpointsinsurface(u,numberpoints,numbersurface)

for k=1:length(numberpoints)
    u = embedpointinsurface(u,numberpoints(k),numbersurface);
end
