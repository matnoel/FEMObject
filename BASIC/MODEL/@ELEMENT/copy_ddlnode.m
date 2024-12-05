function elem=copy_ddlnode(elem,elem2)
% function elem=copy_ddlnode(elem,elem2)


elem.ddlnode = elem2.ddlnode;
elem.ddlgauss = elem2.ddlnode;
elem.ddlnodedual = elem2.ddlnodedual;
elem.ddlgaussdual = elem2.ddlnodedual;
       
elem.nbddlpernode=length(ddlnode);
elem.nbddlpergauss=length(ddlgauss);
elem.nbddl=elem.nbddlpernode*getnbnode(elem);

