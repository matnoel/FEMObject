function M = addelem(M,varargin)
% function M = addelem(M,type,connec,'matnum',matnum,'param',param)
% type   = 'TRI3' ,'BARR', 'BEAM', 'QUA4', ...
% connec :    - connectivite de l'element [noeud1 elem1, noeud2 elem1, ... ; noeud1 elem2, ...]
%        : ou - modele dont on prend la connectivite
%
% matnum : numero du materiau de l'element
% param  : arguments necessaires a la definition des elements
%     param est une propriete commune a tous les elements du groupe ou bien un tableau de cellules de la taille le nombre
%     d'elements ajoutes
%

if nargin==2 & isa(varargin{1},'ELEMENTGEOM')
    elem = varargin{1};
    if getnbelem(elem)>0
    number = M.nbelem + [1:getnbelem(elem)];
    M.nbgroupelem=M.nbgroupelem+1 ;
    M.groupelem{M.nbgroupelem}=elem;
    M.groupelem{M.nbgroupelem}=setnumelem(M.groupelem{M.nbgroupelem},number);
    M.repelemingroupelem=[M.repelemingroupelem;...
      repmat(M.nbgroupelem,length(number),1),[1:length(number)]'];
    M.nbelem = number(end);
    end
    
elseif nargin==2 & isa(varargin{1},'MODEL')
% CAS DU RAJOUT D'UN MODELE    
    M2=varargin{1};
    M = addnode(M,M2);
    for j=1:M2.nbgroupelem
    M = addelem(M,M2.groupelem{j});
    end

elseif nargin>=3 & isa(varargin{2},'MODEL')
   
    M2=varargin{2};
    M = changenodenumber(M,[1:M.nbnode]);    
    M2 = changenodenumber(M2,M.nbnode+[1:M2.nbnode]);
    M = addnode(M,M2);
    
    
    for j=1:M2.nbgroupelem
    connec = getconnec(M2.groupelem{j});
 
    options=cell(1,4);
    options{1}='mat';
    options{2} = getcharin('mat',varargin);
    options{3}='param';
    options{4} = getcharin('param',varargin);

     if isempty(options{2}) & isa(M2.groupelem{j},'ELEMENT')
        options{2} = getmaterial(M2.groupelem{j});
     end
   
    if isempty(options{4}) & isa(M2.groupelem{j},'ELEMENT')
        options{4}  = getparam(M2.groupelem{j});
    end
    
    elemtype = varargin{1} ;

    M=addelem(M,elemtype,connec,options{:});
    
    end
    
    M = unique(M);
    M = changenodenumber(M,[1:M.nbnode]);    

else

    number=M.nbelem+[1:size(varargin{2},1)];   
    if ~isempty(number)
    M.nbgroupelem=M.nbgroupelem+1 ; 

       elemtype = varargin{1} ;
       connec = varargin{2} ;
       mat = getcharin('mat',varargin);
       param = getcharin('param',varargin);
       option = getcharin('option',varargin,' ');
       
      funelem = fcnchk(elemtype);
      M.groupelem{M.nbgroupelem}=...
          funelem(getnode(M.node,unique(connec)),number,connec,mat,option,param);
      M.repelemingroupelem=[M.repelemingroupelem;...
      repmat(M.nbgroupelem,length(number),1),[1:length(number)]'];

        M.nbelem=number(end);
    end
end

%M = concatgroupelem(M);
%M = createddlnode(M);
%M=calc_connec(M);



