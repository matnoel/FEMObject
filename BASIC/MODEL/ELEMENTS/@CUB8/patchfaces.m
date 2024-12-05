function faces = patchfaces(elem,connec)
% function faces = patchfaces(elem,connec)

faces = [1 2 6 5, 2 3 7 6, 3 4 8 7, 4 1 5 8, 1 2 3 4, 5 6 7 8];
if nargin==2
    faces = connec(:,faces)';
    faces = reshape(faces,4,6*size(connec,1))';
else
    faces = reshape(faces,4,6)';
end
