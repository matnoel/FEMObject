function [elemin,elemout,nodeplus,xnodein,xnodeout] = lsdivideelem(elem,ls,node)
lsx = getvalue(ls);
[a,connec] = ismember(getconnec(elem),getnumber(node)) ;
lsx = lsx(connec');
lsx=reshape(lsx,size(connec'));
nbelem = getnbelem(elem);
xi = (lsx(1,:)+lsx(2,:))./ (lsx(1,:)-lsx(2,:)); 
repgin = find(lsx(1,:)<0); % element in a gauche et donc out a droite
repdin = find(lsx(2,:)<0);
repgout = find(lsx(1,:)>=0); % element out a gauche et donc in a droite
repdout = find(lsx(2,:)>=0); 
nodeplus = node ;

xnodein = zerosND(2,1,nbelem);
xnodeout = zerosND(2,1,nbelem);
if ~isempty(repgin)
xnodein(1,1,repgin,1)=-1;
xnodein(2,1,repgin,1)=xi(repgin);
xnodeout(1,1,repgin,1)=xi(repgin);
xnodeout(2,1,repgin,1)=1;
end

if ~isempty(repgout)
xnodein(1,1,repgout,1)=xi(repgout);
xnodein(2,1,repgout,1)=1;
xnodeout(1,1,repgout,1)=-1;
xnodeout(2,1,repgout,1)=xi(repgout);
end

connec=getconnec(elem);
nodeplus = node;
numplus = max(getnumber(node))+[1:nbelem];
xnode = getcoord(node,getconnec(elem)');

xi=MYDOUBLEND(reshape(xi,1,1,nbelem));
x = calc_x(elem,xnode,xi) ; 
x=double(x) ; x = permute(x,[3,2,1]);

nodeplus = addnode(node,x,numplus);
connecin = [connec(repgin,1),numplus(repgin);numplus(repgout),connec(repgout,2)];
connecout = [connec(repgout,1),numplus(repgout);numplus(repgin),connec(repgin,2)];

elemin  = SEG2(nodeplus,getnumber(elem),connecin);
elemout = SEG2(nodeplus,getnumber(elem),connecout);

elemin={elemin};
elemout={elemout};
xnodein = {xnodein};
xnodeout = {xnodeout};



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

end

for k=1:length(elemout)
if isempty(matout)
   elemout{k}=setlstype(elemout{k},'out');
   elemout{k}=setlsnumber(elemout{k},getnumber(ls));
   elemout{k} = setlsenrich(elemout{k},0);   
end
elemout{k} = setmaterial(elemout{k},getmaterial(elem));
end

