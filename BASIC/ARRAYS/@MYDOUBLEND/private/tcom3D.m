function tcomJ = tcom3D(J)
% function tcomJ = tcom3D(J)

dim = size(J,1);
n = size(J,3);
tcomJ = zeros(size(J));
switch dim
    case 1
        tcomJ = ones(1,1,n);
    otherwise
        for i=1:dim
            for j=1:dim
                Jr = J ;
                Jr(:,j,:) = [];
                Jr(i,:,:) = [];
                tcomJ(j,i,:) = (-1)^(i+j)*det3D(Jr);
            end
        end
end
