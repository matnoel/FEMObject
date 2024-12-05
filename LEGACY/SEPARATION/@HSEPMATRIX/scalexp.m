function H=scalexp(H)
% H=scalexp(H)
% Certaines manipulations (evaluation/random) reduisent des dimension
% a un scalaire : il faut faire un peu de menage


H.F=cellfun(@(HF) scalexp(HF),H.F,'uniformoutput',false );
% for d=1:H.dim
%     for r=1:H.m
%         H.F{r,d}=scalexp(H.F{r,d});
%     end
% end

T={};
Hdim=0;
for d=1:H.dim
    if isa(H.F{1,d},'double')
        H.alpha=prod([full(cell2mat(H.F(:,d)))';H.alpha],1);
    elseif isa(H.F{1,d},'SEPMATRIX')
        connect = [0, ones(1,getdim(H.F{1,d})) ];
        T=[T TREE(connect)];
        Hdim=Hdim+1;
        HF(:,Hdim)=H.F(:,d);
    elseif isa(H.F{1,d},'HSEPMATRIX')
        T=[T H.F{1,d}.tree];
        Hdim=Hdim+1;
        HF(:,Hdim)=H.F(:,d);
    end
end

if Hdim>0
    H.dim=Hdim;
    H.tree=TREE(T);
    % Le noeud disparait !
    if H.dim==1
        newH=H.F{1,1}*H.alpha(1);
        for r=2:H.m
            newH=newH+H.F{r,1}*H.alpha(r);
        end
        H=newH;
    end
    
else
    H=sum(H.alpha);
end




