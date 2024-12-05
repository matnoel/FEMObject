function R=setmasse(R,DLmasse)
if ~isempty(DLmasse) && (~isa(DLmasse,'MULTIMATRIX') || length(DLmasse)~=getm(R))
    error('la masse stochastique n''a pas la bonne taille')
end

R.DLmasse=DLmasse;

