function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)
% function [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual] = calc_ddl(mat,elem)

syscoordnode = getsyscoord(elem);
syscoordgauss = getsyscoordlocal(elem);

ddlnode = DDL(DDLVECT('U',syscoordnode,'TRANS'),...
    DDLVECT('R',syscoordnode,'ROTA'));
ddlnodedual = DDL(DDLVECT('F',syscoordnode,'TRANS'),...
    DDLVECT('M',syscoordnode,'ROTA'));

syscoord = getsyscoord(elem);

ddlgauss = DDL(DDLVECT('EPS',syscoordgauss,'TRANS'),...
    DDLVECT('GAM',syscoord,'ROTA'));
ddlgaussdual = DDL(DDLVECT('EFF',syscoordgauss,'TRANS'),...
    DDLVECT('MOM',syscoord,'ROTA'));
