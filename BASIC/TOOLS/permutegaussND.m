function gauss=permutegaussND(gauss)

gauss.coord = MYDOUBLEND(permute(gauss.coord,[4,2,3,1]));
gauss.w     = MYDOUBLEND(permute(gauss.w(:),[2,3,4,1]));

