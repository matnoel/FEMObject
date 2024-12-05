function an=prodscal(a,b,varargin)

if nargin==3 
if isempty(varargin{1})
an = prodscal(a,b);
else
    if ~israndom(a) & ~israndom(b)
       an = a'*expect(varargin{1})*b;
    elseif ~israndom(a)
       a = a'*varargin{1}';
       an = prodscal(a',b); 
    else
       b = varargin{1}*b;
       an = prodscal(a,b);
    end
end
elseif nargin==4

    if ~israndom(a) & ~israndom(b)
       if ~israndom(varargin{1}) 
          an = prodscal(a,b,varargin{1},expect(varargin{2})); 
       elseif ~israndom(varargin{2})
          an = prodscal(a,b,expect(varargin{1}),(varargin{2})); 
       else
        if ~isempty(varargin{1})
        b = varargin{1}*b;
        end
        if ~isempty(varargin{2})
       a = a*varargin{2};
        end   
           an = prodscal(a,b);
       end
       
    elseif ~israndom(a)
       if ~israndom(varargin{1}) & ~israndom(varargin{2})
          an = prodscal(a,expect(b),varargin{1},varargin{2}); 
       else
          an = prodscal(a*varargin{2},varargin{1}*b);
       end
    else
        if ~isempty(varargin{1})
       b = varargin{1}*b;
        end
        if ~isempty(varargin{2})
       a = a*varargin{2};
        end
   
       an = prodscal(a,b);     
    end

else
if ~israndom(a)
    an = a'*expect(b);
elseif ~israndom(b)
    an = expect(a)'*b;
else
    if iscell(a) & ~iscell(b)
    b = mat2cell(b);
    end
    if iscell(b) & ~iscell(a)
    a = mat2cell(a);
    end
    if iscell(a) & iscell(b)
    an = prodscal(getvalue(a.MULTIMATRIX),getvalue(b.MULTIMATRIX));
    else
    an = sum(sum(double(b).*double(a)));    
    end
%    an = sum(sum(double(b).*double(a)));
end

end
    an=full(an);
    

%   an = b.MULTIMATRIX(:)'*a.MULTIMATRIX(:);
%   s = size(an);
%   an = reshape(full(sum(double(an),2)),s); 

