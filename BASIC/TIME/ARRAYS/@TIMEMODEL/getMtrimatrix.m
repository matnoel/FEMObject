function Mt = getMtrimatrix(L)

dt = getdt(L);
nt = getnt(L);
p=getapproxparam(L,'p');

switch getapproxparam(L,'type')
    
    case {'default','CG'}
        Mt=calcul_symbolique(3,p);
        Mt=reshape(Mt(:)*dt,p+1,p+1,p+1,nt);
        [I,J,K,V]=trimatrixelem(Mt,nt,1);
        Mt=SPARSETENSOR([ I J K ], V);
        
    case 'DG'
        matelem = calcul_symbolique(3,p);
        I=[];J=[];K=[];V=[];
        % Peut largement etre optimise....
        for i=1:p+1
        for j=1:p+1
        for k=1:p+1
            I = [I;[i:(p+1):nt*(p+1)]'];
            J = [J;[j:(p+1):nt*(p+1)]'];
            K = [K;[k:(p+1):nt*(p+1)]'];
            V = [V;dt'*matelem(i,j,k)];
        end
        end
        end
        Mt=SPARSETENSOR([ I J K ], V);
        
    otherwise
        error('pas defini')
end




function [i,j,k,val] = trimatrixelem(me,nt,p)
% copie de la fonction trimatrixelem de ELEMENT

numddlelem = [(1:nt)' (2:nt+1)'];
nbddl = p+1;
nbelem = nt ;
s=[nt+1 nt+1 nt+1];

numddlelem=reshape(numddlelem',1,1,nbddl,nbelem);

repcol = repmat(permute(numddlelem,[3,2,1,4]),            [1,nbddl,nbddl,1])+...
    repmat(permute(s(1)*     (numddlelem-1),[1,3,2,4]),   [nbddl,1,nbddl,1])+...
    repmat(        s(1)*s(2)*(numddlelem-1),              [nbddl,nbddl,1,1]);
repcol=  reshape(repcol,nbddl^3,nbelem);


% Assemblage des elements
repu=unique(repcol);
[~,repcolu]=ismember(repcol,repu);
Me=cell(1,nbelem);
Me(:)={spalloc(length(repu),1,nbddl^3)};


for e=1:nbelem
    rep=repcolu(:,e);
    Me{e}(rep)=me(:,:,:,e);
end

[i,j,k]=ind2sub(s,repu);
val=full(sum([Me{:}],2));




