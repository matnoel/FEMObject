function V = shinozuka_cpp(c,z,phi,k,x)
% function V = shinozuka_cpp(c,z,phi,k,x)

[nx,dim] = size(x);
nu = size(k,2);

switch dim
    case 1
        V = shinozuka1D(nu,nx,c,z,phi,k,x);
    case 2
        V = shinozuka2D(nu,nx,c(1,:),c(2,:),z,phi,k(1,:),k(2,:),x(:,1),x(:,2));
    case 3
        V = shinozuka3D(nu,nx,c(1,:),c(2,:),c(3,:),z,phi,k(1,:),k(2,:),k(3,:),x(:,1),x(:,2),x(:,3));
end

end