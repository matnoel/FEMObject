function W = solve_alphaupdate(A,b,W,Wtilde)
%  W = solve_alphaupdate(A,b,W,Wtilde)
%  Optimisation par les alpha_i de W de la fonction coût suivante :
%         J(W.alpha)=|| AW-b ||^2

% Creation des operateurs :
    % WAW_ij = < W_i | W_j >_A
    % WB_i   = < W_i | B >
if nargin>2
    if nargin==3
        Wtilde = W;
    end
    WAW=zeros(W.m);
    WB =zeros(W.m,1);
    for I = 1:Wtilde.m
        for J = 1:W.m
            WAW(I,J) = fastprodscal(truncate(Wtilde,I),truncate(W,J),A) / (Wtilde.alpha(I)*W.alpha(J));
        end
        WB(I) = fastprodscal(truncate(Wtilde,I) , b) / (Wtilde.alpha(I));
    end
    % On resout en alpha :
        % WAW*alpha = WB
    W.alpha = (WAW\WB)';
elseif nargin==2 % Comprendre A=Id
    W=b;
    b=A;
    WAW=zeros(W.m);
    WB =zeros(W.m,1);
    for I = 1:W.m
        for J = 1:W.m
            WAW(I,J) = fastprodscal(truncate(W,I),truncate(W,J)) / (W.alpha(I)*W.alpha(J));
        end
        WB(I) = fastprodscal(truncate(W,I) , b) / (W.alpha(I));
    end
    % On resout en alpha :
        % WAW*alpha = WB
    W.alpha = (WAW\WB)';
end
