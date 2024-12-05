function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)
% function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)

if nargin==1 && isa(mat,'ELEMENT')
    elem = mat;
end

syscoordnode = getsyscoord(elem);
syscoordgauss = getsyscoordlocal(elem);

ddlnode = DDL(DDLVECT('U',syscoordnode));
ddlnodedual = DDL(DDLVECT('F',syscoordnode));

switch getoption(elem)
    case {'BORD','FACE'}
        ddlgauss = DDL(DDLVECT('U',syscoordnode));
        ddlgaussdual = DDL(DDLVECT('F',syscoordnode));
    otherwise
        ddlgauss = DDL(DDLTENS2('EP',syscoordgauss));
        ddlgaussdual = DDL(DDLTENS2('SM',syscoordgauss));
end
