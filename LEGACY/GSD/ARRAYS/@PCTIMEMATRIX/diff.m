function u = diff(u,u0)

if ~isa(u.value,'cell')
if nargin==1
D = getDmatrix(u);
u.value = u.value*D';
elseif nargin==2
D = getDmatrix(u,'ini');
u.value = horzcat(u0,u.value)*D';
    
end
     
else
if nargin==1
D = getDmatrix(u);
u.value = multicellmtimes(u.value,D);    
elseif nargin==2
D = getDmatrix(u,'ini');
u.value = multicellmtimes([{u0} , u.value],D);     
end
end