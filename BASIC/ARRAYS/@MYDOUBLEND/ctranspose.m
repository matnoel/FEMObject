function u = ctranspose(u)
% function u = ctranspose(u)

u.double = conj(permute(u.double,[2,1,3:ndims(u.double)]));
