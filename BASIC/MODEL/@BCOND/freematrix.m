function A = freematrix(BC,A)
% function A = freematrix(BC,A)

if ~isa(BC,'BCOND')
    A = freematrix(A,BC);
else
    if all(gettypes(BC)==0) % dirichlet
        A = A(BC.ddlfree,BC.ddlfree);
    elseif all(gettypes(BC)==1 | gettypes(BC)==0) % periodic
        for i=1:length(BC.BC)
            if BC.BC{i}.type==1
                d2=BC.BC{i}.ddlbloque;
                d1=BC.BC{i}.ddlfree;
                A(:,d1)=A(:,d1)+A(:,d2);
                A(d1,:)=A(d1,:)+A(d2,:);
            end
        end
        A(BC.ddlbloque,:)=[];
        A(:,BC.ddlbloque)=[];
    end
    
end
