function u = transpose(u)
% function u = transpose(u)

u.double = permute(u.double,[2,1,3:ndims(u)]);
