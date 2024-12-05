function u = writefilesol(u,indim,q,file)
% function u = writefilesol(u,indim,q,file)

if nargin==4
    u = setfile(u,file);
end
file = getfilesol(u);

nprec = 2;
q = full(double(q(:)));

fid = fopen(file,'w');
fprintf(fid,'MeshVersionFormatted %u\n\n',nprec);

fprintf(fid,'Dimension %u\n\n',indim);

fprintf(fid,'SolAtVertices\n');
fprintf(fid,'%u\n',size(q,1));
fprintf(fid,'1 1\n\n');

fprintf(fid,'%f \n',q);
fprintf(fid,'\n');
fprintf(fid,'$End');
fclose(fid);
