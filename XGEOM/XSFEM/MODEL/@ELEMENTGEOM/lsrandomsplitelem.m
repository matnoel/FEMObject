function [elemin,elemcut,elemout,repin,repcut,repout,xnodein,xnodecut,xnodeout] = lsrandomsplitelem(elem,ls,node)
% function [elemin,elemcut,elemout,repin,repcut,repout,xnodein,xnodecut,xnodeout] = lsrandomsplitelem(elem,ls,node)

if isa(ls,'LSCRACK')
    warning('utiliser lscracksplitelem')
end

nbs = getfemobjectoptions('nbsamplingsrandomsplit');
RV = RANDVARS(ls);
if nbs>1
%r  = lhsrandom(RV,nbs);
r  = randomwithedges(RV,nbs,1);
else
r = isoprobagrid(RV,nbs);
end
ls = randomeval(ls,r,RV);
ls = lseval(ls,node,'LEVELSET');

[elemin,repin] = lsrandomgetelem(elem,ls,'in',node);
[elemout,repout] = lsrandomgetelem(elem,ls,'out',node);
[elemcut,repcut] = lsrandomgetelem(elem,ls,'cut',node);

if numel(repin)+numel(repout)+numel(repcut)<getnbelem(elem)
    error('pas assez de tirages pour determiner la nature des elements')
end

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



if nargin>2 & nargout>6
xnodein = getcoord(node,getconnec(elemin)');
xnodeout = getcoord(node,getconnec(elemout)');
xnodecut = getcoord(node,getconnec(elemcut)');  
end

