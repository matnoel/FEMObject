function [norms,fieldtype,fieldstorage,fieldddl] = calc_error_deter(group,node,Mref,uref,u1,varargin)
% function [norms,fieldtype,fieldstorage,fieldddl] = calc_error_deter(group,node,Mref,uref,u1,varargin)

fieldtype = 'ddlgaussdual';
fieldstorage = 'center';
fieldddl = getddlgaussdual(group);

norms{1} = zerosND(1,1,getnbelem(group));
norms{2} = zerosND(1,1,getnbelem(group));

for j=1:getnbelem(group)
    % fprintf('elem %d \n',j)
    pourcentage(j,getnbelem(group))
    elem = getelem(group,j);
    n = getcharin('intorder',varargin,1);
    
    connec = getconnec(elem);
    xnode = node(connec);
    xnode = getcoord(xnode,connec');
    ue = localize(elem,u1);
    ueref = localize(elem,uref);
    lspce = Mref.ls(connec);
    [subgaussin,subgaussout] = lssubgauss_oneelem(lspce,n,2);
    if isempty(subgaussin.w)
    else
        detJ = abs(calc_detJ(elem,xnode,subgaussin.coord(1,:)));
        [norms{1}(1,1,j) norms{2}(1,1,j)] = suberrorelem_deter(elem,xnode,subgaussin,ue,ueref,detJ,varargin{:});
    end
end
