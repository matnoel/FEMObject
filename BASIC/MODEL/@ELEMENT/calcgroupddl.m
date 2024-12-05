function groupddl=calcgroupddl(elem)

groupddl.ddlnode = elem.ddlnode ;
groupddl.ddlnodedual = elem.ddlnodedual;
groupddl.node = unique(getconnec(elem));
groupddl.nbddl = elem.nbddlpernode;

