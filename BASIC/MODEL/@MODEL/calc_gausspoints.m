function x = calc_gausspoints(S,choice,varargin)
% function x = calc_gausspoints(S,choice,listenode)

node = getnode(S);
x = cell(getnbgroupelem(S),1);
for i=1:getnbgroupelem(S)
    elem = getgroupelem(S,i);
    xnode = node(elem);
    gauss = calc_gauss(elem,choice);
    xgauss = gauss.coord;
    xi = calc_x(elem,xnode,xgauss);
    x{i} = xi(:,:,:);
end
x = POINT([x{:}]);
