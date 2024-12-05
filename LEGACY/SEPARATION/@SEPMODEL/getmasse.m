function masse=getmasse(SM,dim)
if nargin==1
    masse=SM.masse;
elseif nargin==2
    masse=SM.F{dim}.masse;
end