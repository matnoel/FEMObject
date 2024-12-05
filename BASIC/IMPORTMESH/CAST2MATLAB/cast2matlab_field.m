function [field,node,elem,nb] = cast2matlab_field(file)
% function [field,node,elem,nb] = cast2matlab_field(file)

fid = fopen(file);
if fid==-1
    error(['Le fichier ' file ' n''existe pas']);
end

[node,elem,nb] = cast2matlab_nodeelem(fid);

if nb.chpo>0
    nbsupp=nb.node;
    [field.chpo,field.chponame] = lirechamp(fid,nbsupp,'MAT');
end


if nb.chml>0 
    nbsupp=nb.elem;
    [field.chml,field.chmlname] = lirechamp(fid,nbsupp);
end

fclose(fid);

end
