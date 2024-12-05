function P = calc_P(elem)

Pnode=calc_P(getddlnode(elem),getsyscoord(elem));

s = size2D(Pnode);

P = zerosND([s*getnbnode(elem),sizeND(Pnode)]);

for i=1:getnbnode(elem)
    P((i-1)*s(1)+[1:s(1)],(i-1)*s(2)+[1:s(2)]) = Pnode;
end

