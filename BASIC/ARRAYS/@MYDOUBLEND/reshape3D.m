function u = reshape3D(u)
% function u = reshape3D(u)

u = reshape(u.double,[size2D(u),prod(sizeND(u))]);
