function H=heye(T,s)
%      H=heye(T,s)
% T : arbre de regroupement des variables
% s : ddl par variables
% H : identite sous format hsm
S = SEPMATRIX(length(s),1);
for d=1:length(s)
    S(d)=speye(s(d));
end
H = HSEPMATRIX(S,T);
