function [delta,borne,error] =compute_delta(A,P,U,b)


if nargin == 3
    
    B = A'*(P'*U);
    BTU=B'*U;
    delta=zeros(A.n,1);
    borne=zeros(A.n,1);
    parfor k=1:A.n
        
        Bk=eval(B,k);
        BUk=Bk'*U;
        BBm1=inv(Bk'*Bk);
        BBB = BUk'*BBm1*BUk;
        
        tmp=eig(BBB);
        tmp=min(tmp);
        borne(k)=1/sqrt(tmp);
        delta(k)=sqrt(1-tmp);
        
    end
    
else
    
    
    B = A'*(P'*U);
    BTU=B'*U;
    Pb = P*b;
    delta=zeros(A.n,1);
    borne=zeros(A.n,1);
    error=zeros(A.n,1);
    parfor k=1:A.n
        
        Bk=eval(B,k);
        BUk=Bk'*U;
        BBm1=inv(Bk'*Bk);
        BBB = BUk'*BBm1*BUk;
        
        tmp=eig(BBB);
        tmp=min(tmp);
        borne(k)=1/sqrt(tmp);
        delta(k)=sqrt(1-tmp);
        
        % erreur excacte
        uk=eval_sparse(A,k)\b;
        u_perp= U*(U'*uk);
        ur=U * ( (Bk'*U) \ (U'*eval(Pb,k)) );
        error(k)=norm(uk - ur)/norm(uk-u_perp);
        
    end
    
end

