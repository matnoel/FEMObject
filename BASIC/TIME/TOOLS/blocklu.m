function [L,U] = blocklu(A,Nblock)
% function [L,U] = blocklu(A,Nblock)

Nlocal = size(A,1)/Nblock;
N      = Nlocal;

% matrix to factorize
Aorig = A;
A = cell(Nblock,Nblock);

for i=1:Nblock
    repi = N*(i-1)+(1:N);
    for j=1:Nblock
        repj = N*(j-1)+(1:N);
        A{i,j} = Aorig(repi,repj);
    end
end

% multiplier matrix
L = cell(Nblock,Nblock);

% extract U from the modified A
U= cell(Nblock,Nblock);

% loop over columns
for k=1:Nblock
    
    % check to see if the pivot is zero
    diagA = A{k,k};
    [diagL,diagU] = lu(diagA);
    if(abs(det(diagL))<1e-10)
        disp('warning')
        break; % Gaussian Elimination breaks down
    end
    
    
    L{k,k} = diagL;
    U{k,k} = diagU;
    % compute multipliers
    for i=k+1:Nblock
        L{i,k} = A{i,k}/diagU;
        U{k,i} = diagL\A{k,i};
    end
    
    % apply M to remaining submatrix
    for j=k+1:Nblock
        for i=k+1:Nblock
            A{i,j} = A{i,j} -  L{i,k}*U{k,j};
        end
    end
end



