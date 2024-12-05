function T = assembled_tree(T,celltree)
% T = assembled_tree(T,celltree)
%   -> Assembler les arbres de celltree :
%
%             T
%           / | \
%          /  |  \
%          celltree


nbtrees = length(celltree);
sup_connect=[];
cpt=1;
for i=1:nbtrees
    Ttmp = celltree{i};
    sup_connecti = Ttmp.connect + ones(1,length(Ttmp.connect))*cpt;
    sup_connecti(sup_connecti==cpt)=1;
    cpt=cpt+length(sup_connecti);
    sup_connect = [sup_connect sup_connecti];
end
sup_connect =[0 sup_connect];


% Et voila !
T=TREE(sup_connect);