function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)
% function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)

syscoordnode = getsyscoord(elem);
syscoordgauss = getsyscoordlocal(elem);

ddlnode = DDL(DDLSCAL('U'));
ddlnodedual = DDL(DDLSCAL('QN'));

switch getoption(elem)
    case 'BORD'
        ddlgauss = DDL(DDLSCAL('U'));
        ddlgaussdual = DDL(DDLSCAL('QN'));
    otherwise
        ddlgauss = DDL(DDLVECT('DU',syscoordgauss));
        ddlgaussdual = DDL(DDLVECT('Q',syscoordgauss));
end
