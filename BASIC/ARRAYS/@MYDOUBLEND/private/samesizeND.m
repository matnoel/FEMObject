function [u,v] = samesizeND(u,v)
% function [u,v] = samesizeND(u,v)

u = double(u);
v = double(v);

if ndims(u)>2 || ndims(v)>2
    repunu = find(sizeND(u,1:ndims(v)-2)==1);
    repmatu = ones(1,ndims(v)-2);
    repmatu(repunu) = sizeND(v,repunu);
    u = repmat(u,[1,1,repmatu]);
    
    repunv = find(sizeND(v,1:ndims(u)-2)==1);
    repmatv = ones(1,ndims(u)-2);
    repmatv(repunv) = sizeND(u,repunv);
    v = repmat(v,[1,1,repmatv]);
end