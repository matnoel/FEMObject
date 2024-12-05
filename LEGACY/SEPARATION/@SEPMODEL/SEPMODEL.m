classdef SEPMODEL < SEP
    % L'information du model par dimension est contenue dans une
    % structure de donnee :
    %       -F{d}.model
    %       -F{d}.masse
    %       -F{d}.metric ...
    properties
        mapping={};
        metric =SEP();
        masse  =SEP();
    end
    methods
        function SM = SEPMODEL(M,mapping)
            % S = SEPMODEL(M,varargin)
            % Conteneur de model : permet d'identifier
            % chaque dimension d'un objet SEP
            if nargin==0
                M={};
            elseif nargin==1
                mapping={};
            end
            SM=SM@SEP(M);
            SM.mapping=mapping;
            
            if isempty(SM.F) || ~isa(SM.F{1},'struct')
                SM=init_SEPMODEL(SM);
            end
            
        end
    end
end
