function [xq,yq,z] = computeGrid(x1,x2,fout)
% function [xq,yq,z] = computeGrid(x1,x2,fout)
% See function ksdensity

x = linspace(min(x1),max(x1));
y = linspace(min(x2),max(x2));
[xq,yq] = meshgrid(x,y);
orig_state = warning;
warning('off','all');
z = griddata(x1,x2,fout,xq,yq);
warning(orig_state);

end
