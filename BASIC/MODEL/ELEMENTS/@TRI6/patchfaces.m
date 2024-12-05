function faces = patchfaces(elem,connec)
% function faces = patchfaces(elem,connec)

faces = [1 4 2 5 3 6];
if nargin==2
    faces = connec(:,faces)';
    faces = reshape(faces,6,1*size(connec,1))';
else
    faces = reshape(faces,1,6)';
end
