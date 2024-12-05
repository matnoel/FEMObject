function u = FEELEMFIELD(value,varargin)
% constructeur de la classe FEELEMFIELD
% u = FEELEMFIELD(field,'type',type,'storage',storage,'ddl',ddl) ...
% field : cell array
%    field{i} contient les valeurs du champ
%    sur le groupe d'element i
%    field{i}(:,:,e,l)  contient les valeurs du champ sur le point
%    l de l'element e du groupe d'element i
%    l designe un point de gauss pour un stockage 'gauss', le
%    centre de l'element pour un stockage center, un noeud pour un
%    stockage 'node'
% type : 'ddlgauss' 'ddlgaussdual', ...
% storage : 'node' 'gauss' 'center'
% ddl{i} : objet DDL associe au groupe d'element i
%
% u = FEELEMFIELD(field,MODEL) ...
% dans ce cas field est un double de la taille (M.nbelem * ... * ...  )
% le type de stockage est 'center'

switch class(value)
    case 'FEELEMFIELD'
        u = value;
    case 'cell'
        u.type = getcharin('type',varargin);
        u.storage = getcharin('storage',varargin,'center');
        u.value = value;
        u.ddl = getcharin('ddl',varargin,cell(size(u.value)));
        
        u = class(u,'FEELEMFIELD') ;
    case 'double'
        u.type = getcharin('type',varargin,'scalar');
        u.storage = 'center';
        M = getclassin('MODEL',varargin);
        nbelem = 0;
        s = size(value);
        for i=1:M.nbgroupelem
            rep = nbelem + (1:getnbelem(M.groupelem{i}));
            if size(value,1)==getnbelem(M)
                u.value{i} = reshape(value(rep,:),[length(rep),s(2:end)]);
                u.value{i} = permute(u.value{i},[2,3,1]);
            elseif  size(value,3)==getnbelem(M)
                u.value{i} = value(:,:,rep);
            else
                error('pas prevu')
            end
            
            nbelem = rep(end);
        end
        u.ddl = getcharin('ddl',varargin,cell(size(u.value)));
        
        u = class(u,'FEELEMFIELD') ;
    otherwise
        help FEELEMFIELD
        error('Entrer les bonnes donnees en entree')
end
