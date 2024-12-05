function V = shinozukaRand_cpp(c,z,phi,k,x)
% function V = shinozukaRand_cpp(c,z,phi,k,x)

[nx,dim] = size(x);
order = size(k,2);

switch dim
    case 1
        V = shinozukaRand1D(order,nx,c,z,phi,k,x);
    case 2
        V = shinozukaRand2D(order,nx,c,z,phi,k(1,:),k(2,:),x(:,1),x(:,2));
    case 3
        V = shinozukaRand3D(order,nx,c,z,phi,k(1,:),k(2,:),k(3,:),x(:,1),x(:,2),x(:,3));
end

end