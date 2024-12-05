function u = freevector(BC,u,k)
% function u = freevector(BC,u,k)

if all(gettypes(BC)==0) % dirichlet
    rep = cell(ndims(u),1);
    rep(:)={':'};
    if nargin<3
        k=1;
    end
    rep(k)={BC.ddlfree};
    
    u = u(rep{:});
elseif all(gettypes(BC)==1 | gettypes(BC)==0 ) % periodic
    if nargin==3 && k==2
        u=freevector(BC,u',1)';
        return
    elseif nargin==3 && k~=1
        error('pas programme')
    end
    for i=1:length(BC.BC)
        if BC.BC{i}.type==1
            d2=BC.BC{i}.ddlbloque;
            d1=BC.BC{i}.ddlfree;
            u(d1,:)=u(d1,:)+u(d2,:);
        end
    end
    u(BC.ddlbloque,:)=[];
end


