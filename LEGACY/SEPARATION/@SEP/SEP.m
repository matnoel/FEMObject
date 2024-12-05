classdef SEP
    properties
        dim   = 0;      % Nb de dimensions
        m     = 0;      % Nb de mode
        alpha = [];     % Tableau scalaire 1*dim
        F     = {};     % Tableau cellule  m*dim
    end
    
    methods
        function A = SEP(F,alpha,varargin)
            % function A = SEP(F,alpha)
            % function A = SEP(dim)
            % function A = SEP(HSM)
            if nargin==0
                % Laisser les valeurs d'initialisation
                
            elseif isa(F,'HSEP')
                % Expansion de la HSEP
                A=HSEP2SEP(F);
                
            elseif nargin==1 && (isa(F,'SEP') || isa(F,'struct'))
                % Initialisation a partir d'une (sous)classe SEP,
                % ou d'une struct equivalente
                A.m=F.m;
                A.dim=F.dim;
                A.alpha=F.alpha;
                A.F=F.F;
                
            elseif isa(F,'double')
                % Dimension uniquement
                A.dim = F;
                
            elseif isa(F,'cell')
                % Cellule uniquement : comprendre alpha=1 si non
                % specifie
                A.dim = size(F,2);
                A.m   = size(F,1);
                if nargin==1, A.alpha = ones(1,A.m);
                else          A.alpha = alpha;
                end
                A.F = F;
                
            end
        end
    end
end


