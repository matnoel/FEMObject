function u = setcoord(u,coord)
    coord = MYDOUBLEND(coord);

    if size(coord,1)>1
        coord=coord';
        coord=reshape(coord,[size(coord,1),1,size(coord,2),sizeND(coord)]);
        coord=coord';
    end
    
u.MYDOUBLEND = MYDOUBLEND(coord);

