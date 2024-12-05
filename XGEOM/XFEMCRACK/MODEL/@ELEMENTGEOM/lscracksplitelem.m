function [elemcut,elembicut,elemin,repcut,repbicut,repin,xnodecut,xnodebicut,xnodein] = lscracksplitelem(elem,ls,node)
% function [elemcut,elembicut,elemin,repcut,repbicut,repin,xnodecut,xnodebicut,xnodein] = lscracksplitelem(elem,ls,node)

[elemcut,repcut]= lscrackgetelem(elem,ls,'cut',node);
[elemin,repin]= lscrackgetelem(elem,ls,'in',node);

elemin = setlstype(elemin,'indomain');
elemcut = setlstype(elemcut,'cut');
elemcut = setlsnumber(elemcut,getnumber(ls));
elemcut = setlsenrich(elemcut,isenrichsupport(ls));
elemcut = setlsnature(elemcut,getnature(ls));
if isenrichsupport(ls)
elemcut = setparam(elemcut,'lsenrichtype',getenrichtypesupport(ls));
end
[elembicut,repbicut]= lscrackgetelem(elem,ls,'bicut',node);
if ~isa(elembicut,'cell')
    elembicut = {elembicut};
    repbicut = {repbicut};
end
for k=1:length(elembicut)
elembicut{k} = setlsenrich(elembicut{k},isenrichtip(ls,k));
elembicut{k} = setparam(elembicut{k},'tipnumber',k);
elembicut{k} = setlstype(elembicut{k},'bicut');
elembicut{k} = setlsnumber(elembicut{k},getnumber(ls));  
elembicut{k} = setlsnature(elembicut{k} ,getnature(ls));
if isenrichtip(ls,k)
elembicut{k} = setlsenrichtype(elembicut{k},getenrichtypetip(ls,k));

end
end


if nargout>6
xnodecut = getcoord(node,getconnec(elemcut)');  
xnodein = getcoord(node,getconnec(elemin)');  
xnodebicut = cell(1,length(elembicut));
for k=1:length(elembicut)    
xnodebicut{k} = getcoord(node,getconnec(elembicut{k})');
end
end

