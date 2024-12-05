function Aexp=expand2D(A)
% FONCTION DIFFERENTE DE EXPAND !!!
%      kron(A.V{:,2},A.V{:,1})
% et pas
%      kron(A.V{:,1},A.V{:,2})
m=1;
Aexp = A.alpha(m)*kron(A.F{m,2},A.F{m,1});
for m=2:A.m
    Aexp =Aexp + A.alpha(m)*kron(A.F{m,2},A.F{m,1});
end