function M = addelem(M,A,varargin)
% function M = addelem(M,type,connec,mat,'param',param,'option',option)
% type   = 'TRI3' ,'BARR', 'BEAM', 'QUA4', ...
% connec :    - connectivite de l'element [noeud1 elem1, noeud2 elem1, ... ; noeud1 elem2, ...]
%        : ou - model dont on prend la connectivite
%        : ou - element dont on prend la connectivite
% mat    : materiau de l'element
% param  : argument necessaire a la definition des elements (depend de l'implementation des elements)
% option : voir options possibles pour les differents elements
%
% function M = addelem(M,e,mat,'param',param,'option',option)
% e : model (ou element) dont on prend les elements

% if nargin==2 & isa(A,'ELEMENT')
%     elem = A;
%     if getnbelem(elem)>0
%         M.nbgroupelem=M.nbgroupelem+1 ;
%         M.groupelem{M.nbgroupelem}=elem;
%         M.nbelem = M.nbelem + getnbelem(elem);
%     end
% end

if nargin>=2 && isa(A,'ELEMENT')
    
    elem=A;
    if getnbelem(elem)>0
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        param = getcharin('param',varargin,getparam(elem));
        option = getcharin('option',varargin,getoption(elem));
        elem=setparam(elem,param);
        elem=setoption(elem,option);
        elem=setmaterial(elem,mat);
        M.nbgroupelem=M.nbgroupelem+1 ;
        M.groupelem{M.nbgroupelem}=elem;
        M.nbelem = M.nbelem + getnbelem(elem);
        
        if ~ischarin('norenumelem',varargin)
            M = changeelemnumber(M);
        end
    end
    
elseif isa(A,'MODEL')
    % CAS DU RAJOUT D'UN MODELE
    if ~ischarin('norenum',varargin)
        M = changenodenumber(M,1:M.nbnode);
        A = changenodenumber(A,M.nbnode+(1:A.nbnode));
    end
    M = addnode(M,A);
    
    for j=1:A.nbgroupelem
        M=addelem(M,A.groupelem{j},varargin{:});
    end
    
    M = addfacesofmodel(M,A);
    
    if ~ischarin('duplicate',varargin)
        M = unique(M);
    end
    if ~ischarin('norenum',varargin)
        M = changenodenumber(M,(1:M.nbnode));
    end
    
elseif isa(varargin{1},'MODEL')
    M=addnode(M,varargin{1});
    for i=1:varargin{1}.nbgroupelem
        M=addelem(M,A,varargin{1}.groupelem{i},varargin{2:length(varargin)});
    end
    if ~ischarin('duplicate',varargin)
        M=unique(M);
    end
    
elseif isa(varargin{1},'ELEMENTGEOM')
    elem=varargin{1};
    M=addelem(M,A,getconnec(elem),varargin{2:length(varargin)});
    
else
    funelem = fcnchk(A);
    connec = varargin{1};
    
    number=M.nbelem+[1:size(connec,1)];
    
    if ~isempty(number)
        elem = funelem(getnode(M.node,unique(connec)),number,connec,varargin{2:length(varargin)});
        M = addelem(M,elem) ;
    end
    
end

%M = concatgroupelem(M);
%M = calc_connec(M);
