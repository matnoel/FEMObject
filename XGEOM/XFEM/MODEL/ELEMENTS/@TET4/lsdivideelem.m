function [elemin,elemout,nodeplus,xnodein,xnodeout,nodefield] = lsdivideelem(elem,ls,node,nodefield)

connec = getconnec(elem);
xnode = double(getcoord(node));
nbelem = getnbelem(elem);
[a,conneclocal] = ismember(connec,getnumber(node)) ;
numnode = getnumber(node);
numelem = getnumber(elem);
connectet4in = zeros(0,4);
connectet4out = zeros(0,4);

numelemtet4in=[];
numelemtet4out=[];

xnodetet4in = zeros(0,3);
xnodetet4out = zeros(0,3);


if isa(ls,'LEVELSET')
lsx = getvalue(ls);
elseif isa(ls,'double')
lsx = ls;
else
    error('ls doit etre une levelset ou un double')
end
    
    
for e=1:nbelem
connece = conneclocal(e,:);

[conneclocalin,conneclocalout,xlnodeplus]=lsdivide_oneelem(lsx(connece));%

xnodeplus = calc_x(elem,xnode(connece,:),xlnodeplus);

numplus = max(numnode)+[1:size(xnodeplus,1)];
numnodes = [connece,numplus];
connecin = numnodes(conneclocalin);
connecout = numnodes(conneclocalout);
numnode = [numnode;numplus(:)];
xnode = [xnode;xnodeplus];
xlnodeelem = [nodelocalcoord(elem);xlnodeplus];

if nargout>5 && nargin==4
nodefieldplus = getN(elem,xlnodeplus)*nodefield(connece);
nodefield = [nodefield;zeros(max(numnode)-length(nodefield),1)];
nodefield(numplus) = nodefieldplus;
end

connectet4in = [connectet4in;connecin];
numelemtet4in = [numelemtet4in;repmat(numelem(e),size(connecin,1),1)];
xnodetet4in = [xnodetet4in;xlnodeelem(conneclocalin',:)];


connectet4out = [connectet4out;connecout];
numelemtet4out = [numelemtet4out;repmat(numelem(e),size(connecout,1),1)];
xnodetet4out = [xnodetet4out;xlnodeelem(conneclocalout',:)];


end


nodeplus = NODE(xnode,numnode);

elemin = cell(0,1);
elemout = cell(0,1);
xnodein = cell(0,1);
xnodeout = cell(0,1);


if ~isempty(numelemtet4in)
elemin = [elemin, { TET4(nodeplus,numelemtet4in,connectet4in)}];
xnodetet4in = reshape(xnodetet4in,[size(connectet4in'),3]);
xnodein = [xnodein , {permute(xnodetet4in,[1,3,2])}];
end

if ~isempty(numelemtet4out)
elemout = [elemout, { TET4(nodeplus,numelemtet4out,connectet4out)}];
xnodetet4out = reshape(xnodetet4out,[size(connectet4out'),3]);
xnodeout = [xnodeout , {permute(xnodetet4out,[1,3,2])}];
end

if ~isempty(getmaterial(ls))
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end

for k=1:length(elemin)
   elemin{k}=setlstype(elemin{k},'in');
   elemin{k}=setlsnumber(elemin{k},getnumber(ls));
   elemin{k} = setlsenrich(elemin{k},0);
   elemin{k} = setmaterial(elemin{k},matin);
   xnodein{k} = MYDOUBLEND(xnodein{k});
end

for k=1:length(elemout)
if isempty(matout)
   elemout{k}=setlstype(elemout{k},'out');
   elemout{k}=setlsnumber(elemout{k},getnumber(ls));
   elemout{k} = setlsenrich(elemout{k},0);   
end
elemout{k} = setmaterial(elemout{k},getmaterial(elem));
xnodeout{k} = MYDOUBLEND(xnodeout{k});
end

