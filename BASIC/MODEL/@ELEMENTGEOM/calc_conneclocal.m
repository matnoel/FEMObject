function connec=calc_conneclocal(elem,node)

[a,connec] = ismember(elem.connec,getnumber(node)) ;
