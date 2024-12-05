function v = reshape2D(u)
v=reshape3D(u);
v=permute(v,[2,1,3]);
v=reshape(v,size(v,1),size(v,2)*size(v,3));
v=v';


