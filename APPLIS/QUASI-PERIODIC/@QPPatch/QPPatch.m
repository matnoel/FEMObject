classdef QPPatch
   
    properties( SetAccess = public, GetAccess = public )
        model
        type
        conductivity
        solver
    end
    
    
    methods( Access = public )
        %% Constructor
        
        function patch = QPPatch(varargin)
            p = ImprovedInputParser();
            addParameter(p,'model',[]);
            addParameter(p,'type',[]);
            addParameter(p,'conductivity',[]);
            addParameter(p,'solver',[]);
            parse(p,varargin{:});
            patch = passMatchedArgsToProperties(p,patch);
        end        
        
        %% "Get" methods
        
        % Primary
        
        function model = getModel(P)
            model = P.model ;
        end
                
        function cellNb = getCellNb(P)
            cellNb = getCellNb(getModel(P)) ;
        end
        
        function conductivity = getConductivity(P)
            conductivity = P.conductivity ;
        end
        
        function solver = getSolver(P)
            solver = P.solver ;
        end
        
        function type = getType(p)
            type = p.type;
        end
        
        % Secondary
        
        function flag = isempty(P)
            flag = isempty(getCells(P)) ;
        end
        
        function name = typeName(p)
            switch getType(p)
                case 1
                    name = 'SWIP' ;
                case 2
                    name = 'Neumann-Dirichlet' ;
                case 3
                    name = 'Lagrange' ;
            end
        end
                
        % References
        
        function cellModel = getCellModel(P)
            cellModel = getCellModel(getModel(P)) ;
        end
        
        function cellNum = getCellNum(P)
            cellNum = getCellNum(getModel(P)) ;
        end
        
        function cellNb = getModelCellNb(P)
            cellNb = getCellNb(getModel(P)) ;
        end
           
        %% "Set" methods
        
        function P = setModel(P,model)
           P.model = model ; 
        end
        
        function P = setConductivity(P,conductivity)
            P.conductivity = conductivity ;
        end
        
        function P = setSolver(P,solver)
            P.solver = solver ;
        end
        
        function p = setType(p,newType)
            p.type = newType ;
        end
        
        %% External methods signatures
        
        [u,flux,output] = solve(patch,source,value,flux,tolerance,penalty)
       
    end
    
    methods (Static)
        
        patches = batch(varargin)
        
    end
    
end
