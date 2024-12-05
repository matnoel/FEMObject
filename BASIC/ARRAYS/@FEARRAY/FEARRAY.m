function u = FEARRAY(A,f,d)

% u=FEARRAY(A,f,ddlbloque)  
% A : array 
% f : vecteur indiquant les dimensions Elements finis 
% ddlbloque : ddl libres (tous par defaut)

n=ndims(A);
u.fedims = f;
u.nbddl = size(A,f(1));
if nargin<3
    d=[];
end
u.ddlbloque = d;
u.ddlfree = setdiff([1:u.nbddl]',u.ddlbloque);

u.repbloque = cell(1,n);
u.repbloque(:)={':'};
u.repbloque(f)= {u.ddlbloque} ;
u.repfree = cell(1,n);
u.repfree(:)={':'};
u.repfree(f)= {u.ddlfree} ;

u = class(u,'FEARRAY',MYDOUBLE(A));    
superiorto('MYDOUBLE')
inferiorto('PCARRAY')
inferiorto('PCRADIAL')
