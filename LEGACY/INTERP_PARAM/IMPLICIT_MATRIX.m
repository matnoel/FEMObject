classdef IMPLICIT_MATRIX
    % Class implicit matrix
    
    
    
    properties
        mtimes_right
        mtimes_left
        s
        ctranspose_flag
    end
    
    methods (Access = public)
        
        function x = IMPLICIT_MATRIX(s,mtimes_right,mtimes_left)
            x.s=s;
            x.mtimes_right=mtimes_right;
            x.mtimes_left=mtimes_left;
            x.ctranspose_flag = 0;
        end
        
        
        function y = mtimes(x,y)
            if isa(x,'IMPLICIT_MATRIX') && isa(y,'double')
                if length(y)==1
                    mr = @(t) y*(x*t);
                    ml = @(t) y*(t*x);
                    y  = IMPLICIT_MATRIX(x.s,mr,ml);
                elseif x.ctranspose_flag==0
                    y=x.mtimes_right(y);
                else
                    y=x.mtimes_left(y')';
                end
            elseif isa(x,'double') && isa(y,'IMPLICIT_MATRIX')
                y = (y'*x')';
            else
                error('Not implemented');
            end
        end
        
        function s = size(x)
            s=x.s;
        end
        
        function x = ctranspose(x)
            x.ctranspose_flag=~x.ctranspose_flag;
            x.s=fliplr(x.s);
        end
        
        function z = plus(x,y)
            mr = @(t) x*t + y*t;
            ml = @(t) t*x + t*y;
            z  = IMPLICIT_MATRIX(x.s,mr,ml);
        end
        
        function z = uminus(x)
            mr = @(t) - x*t;
            ml = @(t) - t*x;
            z  = IMPLICIT_MATRIX(x.s,mr,ml);
        end
        
        function z = minus(x,y)
            z  = x + (-y);
        end
        
        function z = mldivide(x,y)
            z = gmres(@(t)x*t,y);
        end
        
    end
    
    
end
