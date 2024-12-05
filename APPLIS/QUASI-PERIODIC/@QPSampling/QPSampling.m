classdef QPSampling
    
    properties
        problem
        generator
        solverChoice
        controlSampler
        controlCoef
        coef = 1
    end
    
    methods
        function s = QPSampling(pb,generator,solverChoice,controlSampler)
            if nargin < 4
                controlSampler = [] ;
                if nargin < 3
                    solverChoice = 2 ;
                end
            end
            s.problem = pb ;
            s.generator = generator ;
            s.solverChoice = solverChoice ;
            s.controlSampler = controlSampler ;
        end
        
        function sample = random(s,nb)
            if isnumeric(nb) && nb < 1
                sample = [] ;
                return
            end
            K = s.generator(nb) ;
            [~,out] = homogenise(s.problem,K,0,0,s.solverChoice) ;
            sample = cellfun(@(c)c([1 end]),out.correctedConductivity,...
                'UniformOutput',false);
            sample = s.coef*cat(1,sample{:}) ;
            if ~isempty(s.controlSampler)
                controlSample = random(s.controlSampler,K) ;
                if isempty(s.controlCoef)
                    c = cov([sample controlSample]) ; 
                    c = c([9 14])./var(controlSample,0,1) ;
                    c = repmat(c,size(sample,1),1);
                else
                    c = s.controlCoef ;
                end
                sample = sample - c.*controlSample ;
            end
        end
        
        function sz = size(s)
            sz = [1 2] ;
        end
    end
    
    methods (Static)
        function s = EimControl(pb,varargin)
            interpolate = @(K) QPSampling.eInterpolate(K,varargin{:}) ;
            s = QPSampling(pb,interpolate,2,[]) ;
        end
        
        function Kx = eInterpolate(K,varargin)
            Kx = cell(size(K)) ;
            for n=1:numel(K)
                [~,Kx{n}] = interpolate(EimInterpolator(double(K{n}),varargin{:})) ;
            end
        end
    end
end