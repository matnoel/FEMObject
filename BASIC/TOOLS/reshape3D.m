function u=reshape3D(u)

u=reshape(u,[size2D(u) , prod(sizeND(u))]);
