function w=vertcat(u,v,varargin)
  
if nargin==1
    w=u;
else
if nargin>500 
  warning('ca ne va pas marcher : le vertcat atteint ses limites de recursivite')  
    if isa(u,'PCMATRIX') & isa(v,'PCMATRIX') &&  all(size(u)==1) & all(size(v)==1)   
[rep,pos] = isclassin('PCMATRIX',varargin);
if length(pos)==length(varargin)
    warning('un peu de bricole')
for i=1:length(varargin)
   rep= rep & all(size(varargin{i})==1);
end
if rep
 w{1} = double(u.MULTIMATRIX);
 w{2} = double(v.MULTIMATRIX);
  for i=1:length(varargin)
    w{i+2} = double(varargin{i}.MULTIMATRIX);  
  end
  w = vertcat(w{:});
  w = PCMATRIX(w,[size(w,1),1],getPC(u));
  return
end
end
    end
end

if isa(u,'PCMATRIX') && isa(v,'PCMATRIX')
 w=u;
 w.MULTIMATRIX = vertcat(u.MULTIMATRIX,v.MULTIMATRIX); 
elseif isa(u,'double')
 try
     n=normest(u);
 catch
     n=norm(u);
 end
 if n==0
    w=v;
    w.MULTIMATRIX = vertcat(u,w.MULTIMATRIX);
 else
    u = expand(u * one(getPC(v)));  
    w = vertcat(u,v);
 end

elseif isa(v,'double')
 try
     n=normest(v);
 catch
     n=norm(v);
 end   
 if n==0
    w=u;
    w.MULTIMATRIX = vertcat(w.MULTIMATRIX,v); 
 else
    v = expand(v * one(getPC(u)));
    w = vertcat(u,v);
 end
else
error('pas défini');    
end

if nargin>2
  w = vertcat(w,varargin{:});
end

end

end