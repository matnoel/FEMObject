function K = FEMATRIX(A,M)

if nargin==1 & isa(A,'FEMATRIX')
    K=A;
elseif nargin==1 & isa(A,'MODEL')
    n = getnbddl(A);
    K = sparse(n,n,10*n);
    K = FEMATRIX(K,A);

elseif nargin==2
    if isa(M,'MODEL') | isa(M,'FEARRAY')
        ddlbloque = getddlbloque(M);
    else
        ddlbloque = M ;   
    end

    A = FEARRAY(A,[1,2],ddlbloque);
    K=struct();
    K=class(K,'FEMATRIX',A);
    superiorto('MYDOUBLE')
    inferiorto('PCARRAY')
end

