classdef affine_matrix
    % Class for implicite matrix of type:
    %
    % A = \sum Phi_i(xi_j) A_i
    %
    % where the LU decomposition of A_i is stored
    
    
    properties
        A        % collection of the A_i
        s        % size
        r        % number of term in the decomposition
        n        % size of the sample
        Phi      % coefficients as function of xi : Phi(i,j)
    end
    
    methods (Access = public)
        
        function x = affine_matrix(A,Phi)
            x.s=size(A{1});
            x.A=A;
            x.Phi=Phi;
            x.r = length(A);
            if length(A)==size(A,2)
                x.A=x.A';
            end
            if size(Phi,1)~=x.r
                x.Phi=x.Phi';
            end
            x.n = size(x.Phi,2);
        end
        
        
        function x = plus(x,y)
            if ~all(x.s==y.s)
                error('Not the same matrix size')
            end
            x.A = [x.A;y.A];
            x.r = x.r + y.r;
            x.Phi = [x.Phi ; y.Phi];
        end
        
        function x = mtimes(x,y)
            if isa(x,'affine_matrix') && isa(y,'double')
                x.A = cellfun( @(A) A*y  ,x.A,'UniformOutput',0);
                
            elseif isa(x,'affine_matrix') && isa(y,'affine_matrix')
                [I,J]=ind2sub([y.r,x.r],1:(x.r*y.r));
                x.A = cellfun( @(X,Y) X*Y  ,x.A(J(:)),y.A(I(:)),'UniformOutput',0);
                x.Phi = y.Phi(I,:).*x.Phi(J,:);
                x.r = x.r*y.r;
                
                
            else
                error('not implemented')
            end
            x.s=size(x.A{1});
        end
        
        function x = uminus(x)
            x.Phi=x.Phi*(-1);
        end
        
        function x = minus(x,y)
            x= x + (-y);
        end
        
        function x=ctranspose(x)
            x.A  = cellfun( @(x) x', x.A,'UniformOutput',0);
            x.s  = size(x.A{1});
        end
        
        function x=normalize_basis(x)
            n = cellfun(@norm,x.A);
            x.A = cellfun(@(x) x/norm(x) ,x.A,'UniformOutput',0);
            x.Phi = spdiags(n,0,length(n),length(n))*x.Phi;
        end
        
        function M = eval(X,xi)
            % evaluation of A(xi) in the point
            
            M = cellfun( @(x) x(:) ,X.A,'UniformOutput',0);
            M = [M{:}]*X.Phi(:,xi);
            M = reshape(M,X.s(1),X.s(2)*length(xi));
            
        end
        
        
        
        function M = eval_sparse(X,xi)

            
%             M = cellfun( @(x) x(:) ,X.A,'UniformOutput',0);
%             M = [M{:}]*sparse(X.Phi(:,xi));
%             M = reshape(M,X.s(1),X.s(2)*length(xi));
            
            a=X.Phi(:,xi);
            M=X.A{1}*a(1);
            for i=2:X.r
                M=M+X.A{i}*a(i);
            end
            
        end
        
        function [n,N]=norm(X,type)
            
            if nargin<2
                type='inf';
            end
            N=zeros(1,size(X.Phi,2));
            parfor i=1:size(X.Phi,2)
                N(i)=norm(eval(X,i));
            end
%             N = cellfun( @(i)norm(eval(X,i)) , num2cell([1:size(X.Phi,2)]') );
            n = norm(N,type);
            
        end
        
        
    end
    
end
