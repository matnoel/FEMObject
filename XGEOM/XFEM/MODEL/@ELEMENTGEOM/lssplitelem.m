function [elemin,elemcut,elemout,repin,repcut,repout,xnodein,xnodecut,xnodeout] = lssplitelem(elem,ls,node)
% function [elemin,elemcut,elemout,repin,repcut,repout,xnodein,xnodecut,xnodeout] = lssplitelem(elem,ls,node)

if isa(ls,'LSCRACK')
    warning('utiliser lscracksplitelem')
end

[elemin,repin] = lsgetelem(elem,ls,'in',node);
[elemout,repout] = lsgetelem(elem,ls,'out',node);
[elemcut,repcut] = lsgetelem(elem,ls,'cut',node);


elemcut = setlstype(elemcut,'cut');
elemcut = setlsnumber(elemcut,getnumber(ls));
elemcut = setlsnature(elemcut,getnature(ls));
elemcut = setlsenrich(elemcut,isenrich(ls));

switch getnature(ls)
    case 'material'
        elemin = setlstype(elemin,'in');
        elemin = setmaterial(elemin,getmaterial(ls));
        elemin = setlsnumber(elemin,getnumber(ls));
        elemin = setlsnature(elemin,getnature(ls));
        
    case 'domain'
        
        elemin = setlstype(elemin,'in');
        elemin = setlsnature(elemin,getnature(ls));
        elemin = setlsnumber(elemin,getnumber(ls));
        elemout = setlstype(elemout,'out');
        elemout = setlsnumber(elemout,getnumber(ls));
        elemout = setlsnature(elemout,getnature(ls));
end



if nargin>2 && nargout>6
    xnodein = getcoord(node,getconnec(elemin)');
    xnodeout = getcoord(node,getconnec(elemout)');
    xnodecut = getcoord(node,getconnec(elemcut)');
end

