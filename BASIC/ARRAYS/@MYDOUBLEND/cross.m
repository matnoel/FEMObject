function w = cross(u,v,k)
% function w = cross(u,v,k)

[u,v] = samesizeND(u,v);
switch k
    case 1
        u = [u;zeros([3-size(u,1),1,sizeND(u)])];
        v = [v;zeros([3-size(v,1),1,sizeND(v)])];
    case 2
        u = [u;zeros([1,3-size(u,2),sizeND(u)])];
        v = [v;zeros([1,3-size(v,2),sizeND(v)])];
end
w = MYDOUBLEND(cross(u,v,k));
