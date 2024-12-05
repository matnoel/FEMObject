function faces = patchfaces(elem,connec)
% function faces = patchfaces(elem,connec)

faces = [1 2 3 5 6 7, 1 2 4 5 9 8, 1 3 4 7 10 8, 2 3 4 6 10 9];
if nargin==2
    faces = connec(:,faces)';
    faces = reshape(faces,6,4*size(connec,1))';
else
    faces = reshape(faces,6,4)';
end
