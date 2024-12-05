function [ddlnode,ddlnodedual,ddlgauss,ddlgaussdual]=getddl(elem,field)
% function [ddlnode,ddlnodedual,ddlgauss,ddlgaussdual] = getddl(elem)
% function ddl = getddl(elem,field)
% field = 'ddlnode', ...

if nargin==2
    ddlnode = eval(['elem.' field]);
else
    ddlnode = elem.ddlnode ;
    ddlnodedual = elem.ddlnodedual ;
    ddlgauss = elem.ddlgauss;
    ddlgaussdual = elem.ddlgaussdual;
end