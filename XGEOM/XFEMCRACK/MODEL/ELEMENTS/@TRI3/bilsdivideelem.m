function [eleminin,eleminout,elemoutin,...
    elemoutout,nodeplus,xnodeinin,xnodeinout,...
    xnodeoutin,xnodeoutout,nodefield] = ...
    bilsdivideelem(elem,ls1,ls2,node,nodefield)

connec = getconnec(elem);
xnode = double(getcoord(node));
nbelem = getnbelem(elem);
[a,conneclocal] = ismember(connec,getnumber(node)) ;
numnode = getnumber(node);
numelem = getnumber(elem);
connectri3inin = zeros(0,3);
connectri3outin = zeros(0,3);
connectri3outout = zeros(0,3);
connectri3inout = zeros(0,3);
numelemtri3inin=[];
numelemtri3outin=[];
numelemtri3inout=[];
numelemtri3outout=[];
xnodetri3inin = zeros(0,2);
xnodetri3outin = zeros(0,2);
xnodetri3inout = zeros(0,2);
xnodetri3outout = zeros(0,2);


if isa(ls1,'LEVELSET')
ls1x = getvalue(ls1);
elseif isa(ls1,'double')
ls1x = ls1;
else
    error('ls doit etre une levelset')
end
if isa(ls2,'LEVELSET')
ls2x = getvalue(ls2);
elseif isa(ls2,'double')
ls2x = ls2;
else
    error('ls doit etre une levelset')
end    
    

for e=1:nbelem
connece = conneclocal(e,:);

[conneclocalinin,conneclocalinout,conneclocaloutin,conneclocaloutout,xlnodeplus]=bilsdivide_oneelem(ls1x(connece),ls2x(connece));%

xnodeplus = calc_x(elem,xnode(connece,:),xlnodeplus);

numplus = max(numnode)+[1:size(xnodeplus,1)];
numnodes = [connece,numplus];
connecinin = numnodes(conneclocalinin);
connecinout = numnodes(conneclocalinout);
connecoutout = numnodes(conneclocaloutout);
connecoutin = numnodes(conneclocaloutin);
numnode = [numnode;numplus(:)];
xnode = [xnode;xnodeplus];
xlnodeelem = [nodelocalcoord(elem);xlnodeplus];

if nargout>5 & nargin==4
nodefieldplus = getN(elem,xlnodeplus)*nodefield(connece);
nodefield = [nodefield;zeros(max(numnode)-length(nodefield),1)];
nodefield(numplus) = nodefieldplus;
end

if size(connecinin,1)>0
connectri3inin = [connectri3inin;connecinin];
numelemtri3inin = [numelemtri3inin;repmat(numelem(e),size(connecinin,1),1)];
xnodetri3inin = [xnodetri3inin;xlnodeelem(conneclocalinin',:)];
end

if size(connecinout,1)>0
connectri3inout = [connectri3inout;connecinout];
numelemtri3inout = [numelemtri3inout;repmat(numelem(e),size(connecinout,1),1)];
xnodetri3inout = [xnodetri3inout;xlnodeelem(conneclocalinout',:)];
end

if size(connecoutin,1)>0
connectri3outin = [connectri3outin;connecoutin];
numelemtri3outin = [numelemtri3outin;repmat(numelem(e),size(connecoutin,1),1)];
xnodetri3outin = [xnodetri3outin;xlnodeelem(conneclocaloutin',:)];
end

if size(connecoutout,1)>0
connectri3outout = [connectri3outout;connecoutout];
numelemtri3outout = [numelemtri3outout;repmat(numelem(e),size(connecoutout,1),1)];
xnodetri3outout = [xnodetri3outout;xlnodeelem(conneclocaloutout',:)];
end





end


nodeplus = NODE(xnode,numnode);

eleminin = cell(0,1);
elemoutin = cell(0,1);
xnodeinin = cell(0,1);
xnodeoutin = cell(0,1);
eleminout = cell(0,1);
elemoutout = cell(0,1);
xnodeinout = cell(0,1);
xnodeoutout = cell(0,1);


if ~isempty(numelemtri3inin)
eleminin = [eleminin, { TRI3(nodeplus,numelemtri3inin,connectri3inin)}];
xnodetri3inin = reshape(xnodetri3inin,[size(connectri3inin'),2]);
xnodeinin = [xnodeinin , {permute(xnodetri3inin,[1,3,2])}];
end
if ~isempty(numelemtri3outout)
elemoutout = [elemoutout, { TRI3(nodeplus,numelemtri3outout,connectri3outout)}];
xnodetri3outout = reshape(xnodetri3outout,[size(connectri3outout'),2]);
xnodeoutout = [xnodeoutout , {permute(xnodetri3outout,[1,3,2])}];
end

if ~isempty(numelemtri3inout)
eleminout = [eleminout, { TRI3(nodeplus,numelemtri3inout,connectri3inout)}];
xnodetri3inout = reshape(xnodetri3inout,[size(connectri3inout'),2]);
xnodeinout = [xnodeinout , {permute(xnodetri3inout,[1,3,2])}];
end

if ~isempty(numelemtri3outin)
elemoutin = [elemoutin, { TRI3(nodeplus,numelemtri3outin,connectri3outin)}];
xnodetri3outin = reshape(xnodetri3outin,[size(connectri3outin'),2]);
xnodeoutin = [xnodeoutin , {permute(xnodetri3outin,[1,3,2])}];
end

