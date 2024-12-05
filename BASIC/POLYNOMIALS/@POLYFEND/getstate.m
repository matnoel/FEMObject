function u=getstate(H,e)

if nargin==1
   e = [1:H.n]; 
end

u=H.state(e);