function M = addnode(M,varargin)
% Ajouter des noeuds a une table de noeuds
% function node = addnode(node,newnode)
% newnode = [1,x1,y1,z1;2,x2,y2,z2;...]  ou  un POINT ou un NODE
% les noeuds sont rajoutes aux preexistants et renumerotes


switch class(varargin{1})
    case 'double'   
table=varargin{1};

if nargin == 3 & getindim(M)==size(table,2)
table = [varargin{2}(:),varargin{1}];
elseif nargin == 3 
table(:,1)=varargin{2};
end

    case 'NODE'
table = gettable(varargin{1});
if nargin == 3
table(:,1)=number;
end
    case {'POINT','MYDOUBLEND'}   
if isa(varargin{1},'POINT')
xcoord= getcoord(varargin{1});
end


xcoord=double(xcoord);
xcoord=permute(xcoord(:,:,:),[1,3,2]);
xcoord=reshape(xcoord,size(xcoord,1)*size(xcoord,2),size(xcoord,3));
if nargin==3
number=varargin{2}(:);
else
number  = max([getnumber(M);0])+[1:size(xcoord,1)]' ;
end
table = [number,xcoord];
    otherwise
help NODE/addnode
error('mauvais argument')
end



%newpoint = POINT(table(:,2:end));
%table=table(:,1:getindim(getsyscoord(newpoint))+1);

table = [gettable(M);table];
[a,b,c]=unique(table(:,1));
table=sortrows(table(b,:),1);

M.number=table(:,1);
M.POINT = POINT(table(:,2:end));
M.nbnode = numel(M.number);

