classdef SPARSETENSOR
    % Equivalent de la fonction SPARSE de matlab, mais pour des ordre >2.
    % Pratique pour stoquer et manipuler les operateur des TRILINFORM par
    % exemple.
    
    properties
        indices= [];
        bornes = [];
        ordre  = 0;
        val    = [];
    end
    methods
        function SP = SPARSETENSOR(indices,val)
            if nargin==2
                SP.indices = indices;
                SP.bornes  = max(SP.indices);
                SP.val     = val;
                SP.ordre   = size(indices,2);
            elseif nargin==1 && isa(indices,'MULTIMATRIX')
                %                 value=getvalue(indices);
                %                 if ~isa(value,'double'), error(''); end
                %                 sm=getsm(indices);
                %                 s=gets(indices);
                %                 if issparse(value)
                %
                %                 else
                %
                %                 end
                %
                %                 salut=1
            end
        end
        function SP = permute(SP,p)
            if length(p)~=SP.ordre
                error('Permutation incorrecte')
            end
            SP.indices=SP.indices(:,p);
            SP.bornes=SP.bornes(:,p);
        end
        function SP = ipermute(SP,p)
            if length(p)~=SP.ordre
                error('Permutation incorrecte')
            end
            ip=p;
            for i=1:SP.ordre
                ip(i)=find(p==i);
            end
            SP.indices=SP.indices(:,ip);
        end
        function SP = flippdim(SP)
            % Effet miroir sur les dimensions : la premiere devient la
            % derniere etc
            p=SP.ordre:-1:1;
            SP=permute(SP,p);
        end
        function SP = SPTENSORtimesVECTOR(SP,v,D)
            % function SP = SPTENSORtimesVECTOR(SP,v,D)
            % Multiplication du SPARSETENSOR SP avec un vecteur v
            % sur la dimension D
            
            
            % Compatibilite des dimensions
            if nargin==2
                D=SP.ordre;
            end
            
            if SP.bornes(D)~=length(v)
                error('Incompatibilite des dimensions')
            end
            SP.val=SP.val.*v(SP.indices(:,D));
            SP.ordre=SP.ordre-1;
            SP.indices(:,D)=[];
            SP.bornes(:,D)=[];
            
            % On est sympa : si ordre 2, on retourne une sparse matrix !
            if SP.ordre==2
                SP=sparse(SP.indices(:,1),SP.indices(:,2),SP.val,SP.bornes(1),SP.bornes(2));
            end
        end
        
    end
    
end
