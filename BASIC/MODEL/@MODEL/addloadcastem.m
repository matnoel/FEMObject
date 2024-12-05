function M =addloadcastem(M,f)

% function M =addloadcastem(M,f)
% f = champ par point (composante equivalente aux ddldual des noeuds)
    
M.load{M.nbload+1}=zeros(M.nbddl,1);

for i=1:M.nbnode
[rep,repchpo]=findddl(M.groupddl{M.repnodeingroupddl(i,1)}.ddlnodedual,f.chponame) ;
M.load{M.nbload+1}(M.groupddl{M.repnodeingroupddl(i,1)}.numddl(M.repnodeingroupddl(i,2),rep))=f.chpo(i,repchpo)    ;
end

M.load{M.nbload+1}=SPACEFUN(M.load{M.nbload+1},'dual');
M.nbload=M.nbload+1;


