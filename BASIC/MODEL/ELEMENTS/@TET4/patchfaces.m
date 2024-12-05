function faces = patchfaces(elem,connec)
% function faces = patchfaces(elem,connec)

faces = [1 2 3, 1 2 4, 1 3 4, 2 3 4];
if nargin==2
    faces = connec(:,faces)';
    faces = reshape(faces,3,4*size(connec,1))';
else
    faces = reshape(3,4)';
end
