function detJ = det3D(J)
% function detJ = det3D(J)

dim = size(J,1);
n = size(J,3);
switch dim
    case 1
        detJ = J;
    case 2
        detJ = calcdet2(J);
    otherwise
        detJ = zeros(1,1,n);
        for i=1:dim
            Jr = J ;
            Jr(:,1,:) = [];
            Jr(i,:,:) = [];
            detJ = detJ + (-1)^(i+1)*det3D(Jr).*J(i,1,:);
        end
end

function detJ = calcdet2(J)
detJ = zeros(1,1,size(J,3));
detJ(1,1,:) = J(1,1,:).*J(2,2,:) - J(2,1,:).*J(1,2,:);
return

