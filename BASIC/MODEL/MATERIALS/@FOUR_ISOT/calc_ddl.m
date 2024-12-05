function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)
% function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)

syscoordnode = getsyscoord(elem);
syscoordgauss = getsyscoordlocal(elem);

ddlnode = DDL(DDLSCAL('T'));
ddlnodedual = DDL(DDLSCAL('QN'));

switch getoption(elem)
    case 'BORD'
        ddlgauss = DDL(DDLSCAL('T'));
        ddlgaussdual = DDL(DDLSCAL('QN'));
    otherwise
        ddlgauss = DDL(DDLVECT('DT',syscoordgauss));
        ddlgaussdual = DDL(DDLVECT('Q',syscoordgauss));
end
