function M = changenodenumber(M,varargin)
% function M = changenodenumber(M,varargin)

if nargin==1
    oldnumber = getnumber(M.node);
    newnumber = (1:M.nbnode)';
elseif nargin==2
    oldnumber = getnumber(M.node);
    newnumber = varargin{1};
elseif nargin==3
    oldnumber = varargin{1};
    newnumber = varargin{2};
end

M.node = changenodenumber(M.node,oldnumber,newnumber);
M.nbnode = getnbnode(M.node) ;

for j=1:M.nbgroupelem
    M.groupelem{j} = changenodenumber(M.groupelem{j},oldnumber,newnumber) ;
end

try
    M = applyfunctiontofaces(M,@changenodenumber,oldnumber,newnumber);
catch
    warning('changenodenumber not applied to faces')
end

%%% cas nargin==3 ????
% currentnumber=getnumber(M.node);
% if length(currentnumber)<oldnumber
%     [rep,loc] = ismember(currentnumber,oldnumber);
%     if length(loc)~=length(oldnumber)
%         error('pas le on nombre de noeuds')
%     end
%     oldnumber = oldnumber(loc);
%     newnumber = newnumber(loc);
%     
% end


