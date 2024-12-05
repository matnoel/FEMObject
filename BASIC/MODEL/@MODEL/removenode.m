function M = removenode(M,numnodeelim,varargin)
% function M = removenode(M,numnodeelim,numnodereplace)
% supprime des noeuds du MODEL M
% numnodeelim : numeros des noeuds a eliminer
% numnodereplace (facultatif) : nouveaux numeros des noeuds elimines pour la connectivite des elements (facultatif)

if nargin==2
    M.node = removenode(M.node,numnodeelim);
    numelem = findelemwithnode(M,numnodeelim);
    if any(numelem)
        disp('suppression des elements contenant les noeuds elimines')
        M=removeelem(M,numelem);
    end
elseif nargin==3
    numnodereplace=varargin{1};
    oldnumber = getnumber(M.node);
    M.node = removenode(M.node,numnodeelim);
    newnumber = oldnumber ;
    [a,b]=ismember(numnodeelim,oldnumber);
    % newnumber(b)=nodereplace;
    % M = changenodenumber(M,oldnumber,newnumber);
    newnumber(b(a))=numnodereplace(a);
    M = changenodenumber(M,oldnumber,newnumber);
end

M.nbnode = getnbnode(M.node);

M = applyfunctiontofaces(M,@removenode,numnodeelim,varargin{:});


