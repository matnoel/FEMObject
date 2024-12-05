function [elemin,elemout,nodeplus,xnodein,xnodeout] = lscrackdivideelem(elem,ls,node)
lssup = getlssupport(ls);
if strcmp(getlstype(elem),'bicut') 
 k = getparam(elem,'tipnumber');
 lstip = getlstip(ls,k);
[elemin,elemout,nodeplus,xnodein,xnodeout,lstipval] = ...
    lsdivideelem(elem,lssup,node,getvalue(lstip));
 lstip = setvalue(lstip,lstipval);
else
[elemin,elemout,nodeplus,xnodein,xnodeout] = ...
    lsdivideelem(elem,lssup,node);   
end

elemin = [elemin,elemout];
xnodein = [xnodein,xnodeout];
elemout = {};
xnodeout={};

if strcmp(getlstype(elem),'bicut') 
 
elemtemp = {};
xnodetemp = {};

for l=1:length(elemin)    
elemtemp = [elemtemp,{lsgetelem(elemin{l},lstip,'in',nodeplus)}];
elemtemp = [elemtemp,{lsgetelem(elemin{l},lstip,'out',nodeplus)}];
[elemsubin,elemsubout,nodeplus,xnodesubin,xnodesubout] = ...
    lsdivideelem(lsgetelem(elemin{l},lstip,'cut',nodeplus),lstip,nodeplus); 
elemtemp = [elemtemp,elemsubin,elemsubout];
xnodetemp = [xnodetemp,xnodesubin,xnodesubout];
end

elemin = elemtemp;
xnodein = xnodetemp;
end

for k=1:length(elemin)
elemin{k} = setlstype(elemin{k},'indomain');
end


%for k=1:length(elemin)
%calc_N(elemin{k},'nbddlpernode',1)

