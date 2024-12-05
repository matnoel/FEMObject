function M = assemble_trimatrixelem(S,me,varargin)
% s taille de la matrice finale
type='SPARSETENSOR';
% if (~isa(varargin{1},'char'))
%     % Par defaut
%     type='MULTIMATRIX';
% elseif strcmp(varargin{1},'MULTIMATRIX')||...
%        strcmp(varargin{1},'CELLSPARSE')||...
%        strcmp(varargin{1},'SPARSETENSOR')
%     type=varargin{1};
%     varargin=varargin(2:end);
% end

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
s=[S.nbddl,S.nbddl,S.nbddl];
i=[];j=[];k=[];val = [];
for p=liste
    [ip,jp,kp,valp]=trimatrixelem(S.groupelem{p},double(me{p}),s);
    i =[i;ip];
    j =[j;jp];
    k =[k;kp];
    val=[val;full(valp)];
end
    
% FORMAT DE SORTIE
if strcmp(type,'MULTIMATRIX')
    % MULTIMATRIX VERSION
    
    ki= (k-1)*S.nbddl+i ;
    M=MULTIMATRIX(sparse(ki,j,val,s(1)^2,s(2)),[S.nbddl S.nbddl],[S.nbddl 1]);
    
elseif strcmp(type,'CELLSPARSE')
    % CELLSPARSE VERSION
    
    M=repmat({sparse([],[],[],S.nbddl,S.nbddl)},S.nbddl,1);
    for ii=1:length(i)
        M{i(ii)}(j(ii),k(ii))=val(ii);
    end
    % II=num2cell(i);
    % JJ=num2cell(j);
    % KK=num2cell(k);
    % VVAL=num2cell(val);
    % cellfun(@(I,J,K,VAL) subsasgn(M,substruct('{}',{I},'()',{J,K}), VAL
    % ),II(:),JJ(:),KK(:),VVAL(:));
    
elseif strcmp(type,'SPARSETENSOR')
    % SPARSETENSOR VERSION
    M=SPARSETENSOR([i j k],val);
end
