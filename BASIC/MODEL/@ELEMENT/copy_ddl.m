function elem=copy_ddl(elem,elem2)
% function elem=copy_ddl(elem,elem2)

elem.ddlnode = elem2.ddlnode;
elem.ddlnodedual = elem2.ddlnodedual;

if strcmp(getoption(elem),'BORD') || strcmp(getoption(elem),'FACE')
elem.ddlgauss = elem2.ddlgauss;
elem.ddlgaussdual = elem2.ddlgaussdual;
else
elem.ddlgauss = elem2.ddlnode;
elem.ddlgaussdual = elem2.ddlnodedual;    
end

elem.nbddlpernode=length(elem.ddlnode);
elem.nbddlpergauss=length(elem.ddlgauss);
elem.nbddl=elem.nbddlpernode*getnbnode(elem);
