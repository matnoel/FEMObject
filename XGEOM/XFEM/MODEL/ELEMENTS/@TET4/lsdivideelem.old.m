function [elemin,elemout,nodeplus,xnodein,xnodeout,nodefield] = ...
    lsdivideelem(elem,ls,node,nodefield)

conneclocal = calc_conneclocal(elem,node);
xnode = double(getcoord(node));
nbelem = getnbelem(elem);
numnode = getnumber(node);
numelem = getnumber(elem);
connectri3in = zeros(0,3);
connectri3out = zeros(0,3);
connecqua4in = zeros(0,4);
connecqua4out = zeros(0,4);
numelemqua4in=[];
numelemqua4out=[];
numelemtri3in=[];
numelemtri3out=[];

xnodetri3in = zeros(0,2);
xnodequa4in = zeros(0,2);
xnodetri3out = zeros(0,2);
xnodequa4out = zeros(0,2);


if isa(ls,'LEVELSET')
lsx = getvalue(ls);
elseif isa(ls,'double')
%warning('si on rentre un double et pas une levelset, le matreiau n''est pas attribue aux elements')
error('rentrer une LEVELSET')
lsx = ls;
else
    error('ls doit etre une levelset')
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

if nargout>5 & nargin==4
nodefieldplus = getN(elem,xlnodeplus)*nodefield(connece);
nodefield = [nodefield;zeros(max(numnode)-length(nodefield),1)];
nodefield(numplus) = nodefieldplus;
end

switch size(connecin,2)
    case 3
connectri3in = [connectri3in;connecin];
numelemtri3in = [numelemtri3in;repmat(numelem(e),size(connecin,1),1)];
xnodetri3in = [xnodetri3in;xlnodeelem(conneclocalin',:)];

    case 4
connecqua4in = [connecqua4in;connecin];
numelemqua4in = [numelemqua4in;repmat(numelem(e),size(connecin,1),1)];
xnodequa4in = [xnodequa4in;xlnodeelem(conneclocalin',:)];
end

switch size(connecout,2)
    case 3
connectri3out = [connectri3out;connecout];
numelemtri3out = [numelemtri3out;repmat(numelem(e),size(connecout,1),1)];
xnodetri3out = [xnodetri3out;xlnodeelem(conneclocalout',:)];
    case 4
connecqua4out = [connecqua4out;connecout];
numelemqua4out = [numelemqua4out;repmat(numelem(e),size(connecout,1),1)];
xnodeinqua4out = [xnodeinqua4out;xlnodeelem(conneclocalout',:)];
end

end


nodeplus = NODE(xnode,numnode);

elemin = cell(0,1);
elemout = cell(0,1);
xnodein = cell(0,1);
xnodeout = cell(0,1);


if ~isempty(numelemtri3in)
elemin = [elemin, { TRI3(nodeplus,numelemtri3in,connectri3in)}];
xnodetri3in = reshape(xnodetri3in,[size(connectri3in'),2]);
xnodein = [xnodein , {permute(xnodetri3in,[1,3,2])}];
end
if ~isempty(numelemqua4in)
elemin = [elemin, { QUA4(nodeplus,numelemqua4in,connecqua4in)}];
xnodequa4in = reshape(xnodequa4in,[size(connecqua4in'),2]);
xnodein = [xnodein , {permute(xnodequa4in,[1,3,2])}];
end
if ~isempty(numelemtri3out)
elemout = [elemout, { TRI3(nodeplus,numelemtri3out,connectri3out)}];
xnodetri3out = reshape(xnodetri3out,[size(connectri3out'),2]);
xnodeout = [xnodeout , {permute(xnodetri3out,[1,3,2])}];
end
if ~isempty(numelemqua4out)
elemout = [elemout, { QUA4(nodeplus,numelemqua4out,connecqua4out)}];
xnodequa4out = reshape(xnodequa4out,[size(connecqua4out'),2]);
xnodeout = [xnodeout , {permute(xnodequa4out,[1,3,2])}];
end

if isa(ls,'LEVELSET')
switch getnature(ls)
case 'material'
    matin = getmaterial(ls);
    matout = getmaterial(elem);
case 'domain'
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end
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

