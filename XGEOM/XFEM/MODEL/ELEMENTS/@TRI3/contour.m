function [coseg1,node1,contourfield]=contour(elem,node,nodeval,contourval,nodefield)
tolls = getfemobjectoptions('tolerancelevelset');
xnode=double(getcoord(node));
dim=size(xnode,2);
connec = calc_conneclocal(elem,node);


dim=size(xnode,2);
con=nodeval(connec)-contourval;
con=reshape(con,size(connec));
repzero=find(abs(con)<tolls);
con(repzero)=0;

seg=zeros(0,2);
for k=1:3
seg=[seg;connec(:,[k,mod(k,3)+1])];
end

segcut=zeros(size(connec,1),3);
for i=1:3
segcut(:,i) = (con(:,i).*con(:,mod(i,3)+1));   
end
%segcut(find(abs(segcut)<tolls))=0;
segcut=sign(segcut);

rep = find(sum(segcut==0 | segcut==-1,2)>1);
segcut = segcut(rep,:);
concut = sign(con(rep,:));

elemin = find(sum(sign(con)==-1  | sign(con)==0,2)==3);

node1=zeros(0,dim);
nbnode=size(node1,1);
coseg1 = zeros(0,2);

% GESTION DES ELEMENTS DE TYPE I
% la level set coupe deux segments 
elem1 = find((sum(segcut==-1,2)==2) );
if ~isempty(elem1)
segcut1 = segcut(elem1,:);
repseg1=find((segcut1==-1)');
repseg1=repseg1(:);
repseg1=reshape(mod(repseg1(:)-1,3)+1,2,length(elem1))';

repseg11=((repseg1(:,1)-1)*size(connec,1)+rep(elem1));
repseg12=mod(repseg11+size(connec,1)-1,numel(connec))+1;
xi1 = con(repseg11)./(con(repseg11)-con(repseg12));
if nargout==3
valP1 = (1-xi1).*nodefield(connec(repseg11),:)+xi1.*nodefield(connec(repseg12),:);
end
xi1=repmat(xi1,1,dim);
P1 = (1-xi1).*xnode(connec(repseg11),:)+xi1.*xnode(connec(repseg12),:);
repseg11=((repseg1(:,2)-1)*size(connec,1)+rep(elem1));
repseg12=mod(repseg11+size(connec,1)-1,numel(connec))+1;
xi1 = con(repseg11)./(con(repseg11)-con(repseg12));
if nargout==3
valP2 = (1-xi1).*nodefield(connec(repseg11),:)+xi1.*nodefield(connec(repseg12),:);
end
xi1=repmat(xi1,1,dim);
P2 = (1-xi1).*xnode(connec(repseg11),:)+xi1.*xnode(connec(repseg12),:);

node1 = [node1;P1;P2];
coseg1 = [coseg1;...
 [[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);
if nargout==3
contourfield = [contourfield;valP1;valP2];
end
end



% GESTION DES ELEMENTS DE TYPE II
% la level set coupe un segment et passe par un point
elem1 = find((sum(segcut==-1,2)==1) & (sum(concut==0,2)==1));
if ~isempty(elem1)
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');
reppo1 = reppo1(:);
reppo1 = mod(reppo1(:)-1,3)+1;
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
if nargout==3
valP1 = nodefield(connec(reppo11),:);
end
segcut1 = segcut(elem1,:);
repseg1=find((segcut1==-1)');
repseg1=repseg1(:);
repseg1=mod(repseg1(:)-1,3)+1;
repseg11=((repseg1(:,1)-1)*size(connec,1)+rep(elem1));
repseg12=mod(repseg11+size(connec,1)-1,numel(connec))+1;
xi1 = con(repseg11)./(con(repseg11)-con(repseg12));
if nargout==3
valP2 = (1-xi1).*nodefield(connec(repseg11),:)+xi1.*nodefield(connec(repseg12),:);
end
xi1=repmat(xi1,1,dim);
P2 = (1-xi1).*xnode(connec(repseg11),:)+xi1.*xnode(connec(repseg12),:);

node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);
if nargout==3
contourfield = [contourfield;valP1;valP2];
end
end

% GESTION DES ELEMENTS DE TYPE III 
% la level passe par 2 noeuds opposes
elem1 = find((sum(segcut==0,2)==3) & (sum(concut==0,2)==2));
if ~isempty(elem1)
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');

reppo1=reshape(mod(reppo1(:)-1,3)+1,2,length(elem1))';
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
reppo12 = ((reppo1(:,2)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
P2=xnode(connec(reppo12),:);
if nargout==3
valP1 = nodefield(connec(reppo11),:);
valP2 = nodefield(connec(reppo12),:);
end
node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);
if nargout==3
contourfield = [contourfield;valP1;valP2];
end
end



% GESTION DES ELEMENTS DE TYPE IV
% la level passe par 2 noeuds contigus
elem1 = find((sum(segcut==0,2)==3) & (sum(concut==0,2)==2));
if ~isempty(elem1)
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');
reppo1=reshape(mod(reppo1(:)-1,3)+1,2,length(elem1))';
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
reppo12 = ((reppo1(:,2)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
P2=xnode(connec(reppo12),:);
if nargout==3
valP1 = nodefield(connec(reppo11),:);
valP2 = nodefield(connec(reppo12),:);
end

node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);
if nargout==3
contourfield = [contourfield;valP1;valP2];
end
end

% GESTION DES ELEMENTS DE TYPE V
% la level passe par 3 noeuds contigus
elem1 = find((sum(concut==0,2)==3));
if ~isempty(elem1)
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');
reppo1=reshape(mod(reppo1(:)-1,3)+1,3,length(elem1))';
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
reppo12 = ((reppo1(:,2)-1)*size(connec,1)+rep(elem1));
reppo13 = ((reppo1(:,3)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
P2=xnode(connec(reppo12),:);
P3=xnode(connec(reppo13),:);
if nargout==3
valP1 = nodefield(connec(reppo11),:);
valP2 = nodefield(connec(reppo12),:);
valP3 = nodefield(connec(reppo13),:);
end

node1 = [node1;P1;P2;P3];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]';...
                   [length(elem1)+1:2*length(elem1)]',[2*length(elem1)+1:3*length(elem1)]']];
nbnode=size(node1,1);
if nargout==3
contourfield = [contourfield;valP1;valP2;valP3];
end
end


% GESTION DES ELEMENTS DE TYPE V
% la level fait le tour de l'element
elem1 = find((sum(concut==0,2)==3));
if ~isempty(elem1)
for i=1:3
P1 = xnode(connec(rep(elem1),i),:);
P2 = xnode(connec(rep(elem1),mod(i,3)+1),:);
if nargout==3
valP1 = nodefield(connec(reppo11),:);
valP2 = nodefield(connec(reppo12),:);
end

node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);
if nargout==3
contourfield = [contourfield;valP1;valP2];
end

end    
end



node1=NODE(node1,1:size(node1,1));
coseg1 = SEG2(node1,1:size(coseg1,1),coseg1);

