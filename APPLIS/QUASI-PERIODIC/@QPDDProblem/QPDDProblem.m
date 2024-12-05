classdef QPDDProblem < QPProblem
    
   properties
       patch % QPProblem array
       patchCells % cell array of double arrays
       maxIterationsGL % integer. Maximum iterations for global-local algorithm
       verboseGL % boolean. Verbose for global-local algorithm
       coordTol % scalar double
   end
   
   methods (Access = public)
       %% Constructor methods
       function pb = QPDDProblem(model,conductivity,source,bc,patch,...
               patchCells,varargin)
           p = ImprovedInputParser();
           addParameter(p,'maxIterationsGL',20,@isscalar);
           addParameter(p,'verboseGL',true,@islogical);
           addParameter(p,'coordTol',1e-9,@isscalar);
           parse(p,varargin{:});
           pb = pb@QPProblem(model,conductivity,source,bc,p.unmatchedArgs{:});
           pb = passMatchedArgsToProperties(p,pb);
           pb.patch = patch ;
           pb.patchCells = patchCells ;
       end
       
       %% Sundry
       
       function cellList = allPatchCells(pb)
           % cellList = allPatchCells(pb)
          cellList = pb.patchCells ;
          cellList = cat(1,cellList{:}) ;
       end
       
       %% External methods signatures
       
       complete = assembleGlobalLocal(pb,glob,loc,fullModel)
       
       [feModel,mesherTime] = buildFEModel(pb,method)
       
       K = getCompleteK(pb,fullModel)
       
       source = getCompleteSource(pb,fullModel)
       
       [sol,output] = solve(pb,localSolvers,tolerance,reference,relaxation,fullModel)
   end
   
   methods (Static)
       
   end
end