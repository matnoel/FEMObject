classdef implicit_affine_matrix
    % Class for implicite matrix of type:
    %
    % A = \sum Phi_i(xi_j) A_i^{-1}
    %
    % where the LU decomposition of A_i is stored
    
    
    properties
        L        % left part
        U        % right part
        s        % size
        r        % number of term in the decomposition
        n        % size of the sample
        Phi      % coefficients as function of xi : Phi(i,j)
    end
    
    methods (Access = public)
        
        function x = implicit_affine_matrix(A,Phi)
            x.s=size(A{1});
            x.Phi=Phi;
            x.r = length(A);
            if size(Phi,1)~=x.r
                x.Phi=x.Phi';
            end
            x.L = cell(x.r,1);
            x.U = cell(x.r,1);
            for i=1:x.r
                [x.L{i},x.U{i}]=lu(A{i});
            end
            x.n=size(x.Phi,2);
        end
        
        
        function x = plus(x,y)
            if prod(x.s~=y.s)
                error('Not the same matrix size')
            end
            x.L = [x.L;y.L];
            x.U = [x.U;y.U];
            x.r = x.r + y.r;
            x.Phi = [x.Phi , y.Phi];
        end
        
        
        function x = mtimes(x,y)
            if isa(x,'implicit_affine_matrix') && isa(y,'double')
                
                ULy = cellfun( @(L,U) U\(L\y) ,x.L,x.U,'UniformOutput',0);
                x = affine_matrix( ULy,x.Phi );
                
            elseif isa(x,'implicit_affine_matrix') && isa(y,'affine_matrix')
                
                x = mtimes_truncation(x,y);
                
            end
        end
        
        function x = affine_matrix(x)
            x=x*speye(x.s);
        end
        
        function x = uminus(x)
            x.Phi=x.Phi*(-1);
        end
        
        function x = minus(x,y)
            x= x + (-y);
        end
        
        function x=ctranspose(x)
            Ltmp = cellfun( @(x) x', x.U,'UniformOutput',0);
            x.U  = cellfun( @(x) x', x.L,'UniformOutput',0);
            x.L=Ltmp;
        end
        
        
%         function x = subsref(x,S)
%             
%             
%             if length(S.subs{1})== x.r
%                 x.a = S.subs{1};
%             else
%                 error('Not the good length')
%             end
%             
%         end
        
    end
    
    
end
