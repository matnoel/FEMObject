function A=fundim(A,fun,dim)
% function A=fundim(A,fun,dim)
% Application de la fonction 'fun' sur l'ensemble des
% dimensions dim
A.F(:,dim)=cellfun(fun,A.F(:,dim),'UniformOutput',0);