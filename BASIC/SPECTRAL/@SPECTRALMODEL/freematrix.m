function A = freematrix(S,A)
% function A = freematrix(S,A)

if isa(A,'SPECTRALMODEL')
    A = freematrix(A,S);
else
    if ~isempty(S.bc)
        switch S.bc.type
            case 'dirichlet'
                rep = S.bc.node;
                A(:,rep)=[];
                A(rep,:)=[];
            case 'periodic'
                if S.dim==1
                    A(:,1) = A(:,1)+A(:,end);
                    A(1,:) = A(1,:) + A(end,:);
                    A(:,end)=[];
                    A(end,:)=[];
                else
                    
                end
        end
    end
end