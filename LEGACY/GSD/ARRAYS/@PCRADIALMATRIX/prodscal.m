function an=prodscal(a,b,varargin)

%an = prodscal(expand(a),expand(b));
if nargin==3  
if isempty(varargin{1})
an = prodscal(a,b);
else
    if ~israndom(a) & ~israndom(b)
       an = prodscal(a,b,expect(varargin{1}));
    elseif ~israndom(a)
       if ~israndom(varargin{1}) 
       an = prodscal(a,expect(b),varargin{1});     
       else
       a = a'*varargin{1}';
       an = prodscal(a',b); 
       end
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
       if ~isempty(varargin{1}) 
       b = varargin{1}*b;
       end
       if ~isempty(varargin{2})
       a = a*varargin{2};
       end
          an = prodscal(a,b);
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
elseif isa(a,'PCMATRIX')
    an = sum(sum(expectmtimes(a,b)));
    %an = prodscal(a,expand(b));
elseif isa(b,'PCMATRIX')
    an = sum(sum(expectmtimes(a,b)));
    %an = prodscal(expand(a),b);    
elseif isa(a,'PCRADIALMATRIX') & isa(b,'PCRADIALMATRIX')
    vn = double(prodscal(a.V,multitranspose(b.V)));
    vn = reshape(vn,a.m,b.m);
    
    ln = expectmtimes(a.L,b.L');
    
    an = sum(sum(ln.*vn));
    
else
    error('pas programme')
end

end
    
    an=full(an);
