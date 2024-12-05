function u = appendnodedata(u,q,metrictype,file)
% function u = appendnodedata(u,q,metrictype,file)

if nargin<3
    metrictype = 1;
end

if nargin==4
    u = setfile(u,file);
end
file = getfilemsh(u);

[~,name,ext] = fileparts(file);

q = full(double(q(:)));
nbnodes = size(q,1);

fid = fopen(file,'a');
fprintf(fid,'$NodeData\n');
fprintf(fid,'1\n');
fprintf(fid,'%s%s:metric\n',name,ext);
fprintf(fid,'1\n');
fprintf(fid,'0.0\n');
fprintf(fid,'3\n');
fprintf(fid,'0\n');
fprintf(fid,'%u\n',metrictype);
fprintf(fid,'%u\n',nbnodes);
fprintf(fid,'%u %f\n',[(1:nbnodes)',q]');
fprintf(fid,'$EndNodeData');
fclose(fid);
