classdef SEPMATRIX < SEP
    methods
        function A = SEPMATRIX(F,alpha,varargin)
            % function A = SEPMATRIX(F,alpha)
            % function A = SEPMATRIX(dim,m)
            % function A = SEPMATRIX(HSM)
            A=A@SEP();
            if nargin==0
                A.dim = 0;
                A.m = 0;
                A.alpha = [];
                A.F = {};
                
            elseif nargin==1 && (isa(F,'struct')||isa(F,'SEP')) && ~isa(F,'HSEP')
                A.dim = F.dim;
                A.m = F.m;
                A.alpha = F.alpha;
                A.F = F.F;
                
            elseif isa(F,'double')
                
                A.dim = F;
                if nargin==1
                    alpha = [];
                end
                A.m = length(alpha);
                A.alpha = alpha;
                A.F = cell(A.m,A.dim);
                
            elseif isa(F,'cell')
                A.dim = size(F,2);
                A.m = size(F,1);
                if nargin==1
                    A.alpha = ones(1,A.m);
                else
                    A.alpha = alpha(:)';
                end
                A.F = F;
                
            elseif isa(F,'TSEPMATRIX')
                A=SEPMATRIX();
                A.dim = F.dim;
                A.m = prod(F.m);
                A.F = cell(A.m,A.dim);
                A.alpha = zeros(1,A.m);
                c = 1;
                alpha = double(F.alpha);
                A.alpha = alpha(:)';
                switch F.dim
                    case 2
                        for j=1:F.m(2)
                            for i=1:F.m(1)
                                A.F{c,1}=F.F{1}{i};
                                A.F{c,2}=F.F{2}{j};
                                c=c+1;
                            end
                        end
                    case 3
                        for k=1:F.m(3)
                            for j=1:F.m(2)
                                for i=1:F.m(1)
                                    A.F{c,1}=F.F{1}{i};
                                    A.F{c,2}=F.F{2}{j};
                                    A.F{c,3}=F.F{3}{k};
                                    c=c+1;
                                end
                            end
                        end
                    case 4
                        for l=1:F.m(4)
                            for k=1:F.m(3)
                                for j=1:F.m(2)
                                    for i=1:F.m(1)
                                        A.F{c,1}=F.F{1}{i};
                                        A.F{c,2}=F.F{2}{j};
                                        A.F{c,3}=F.F{3}{k};
                                        A.F{c,4}=F.F{4}{l};
                                        c=c+1;
                                    end
                                end
                            end
                        end
                    otherwise
                        error('not implemented');
                end
                
            elseif isa(F,'ktensor')
                A.dim = ndims(F);
                A.m = numel(F.lambda);
                A.alpha = F.lambda';
                A.F = cell(A.m,A.dim);
                for j=1:A.dim
                    for i=1:A.m
                        A.F{i,j}=F.U{j}(:,i);
                    end
                end
                
            elseif isa(F,'HSEPMATRIX')
                A = SEPMATRIX(HSEP2SEP(F));
            end
            
            
            % superiorto('SEPSOLVER')
            
        end
    end
end
