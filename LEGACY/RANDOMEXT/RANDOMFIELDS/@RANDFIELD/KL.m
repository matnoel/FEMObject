function RFpc = KL(RF,varargin)
% function RFpc = KL(RF,m,S)
% RF : RANDFIELD
% S : MODEL de type MAILLAGE
% m : ordre de la decomposition
% RFpc : DISCRANDFIELD
%
% function RFpc = KL(RF,m,S,'order',p)
% p : ordre du chaos polynomial pour la decomposition des variables
% aleatoires
%
disp('---- Decomposition de Karhunen-Loeve ----')
S=getclassin('MODEL',varargin);

S=createddlnode(S,DDL('U'));

D=calc_massgeom(S);
x = getcoord(S.node);

if RF.selfcorrel
    disp('Compute correlation matrix ...')
    R=correl(RF,x,x);
    s=std(RF,x);
    R = R.*(s*s');
    opt.disp=0;
    disp('Solve eigenvalue problem ...')
    m=varargin{1};
    [V,L]=eigs(R*D,m+1,'LM',opt);
    V = V/diag(sqrt(sum(V.*(D*V),1)));
    lam=diag(L);
    [lam,ia]=sort(lam,'descend');
    error = sqrt(lam(end)/sum(lam));
    lam = lam(1:m);
    V=V(:,ia(1:m));
    
    disp(['KL truncation for the gaussian field : order ' num2str(m) ])
    fprintf('error estimate for the truncation = %3d\n',error);
    
    if ischarin('rescale',varargin)
        V = V.*repmat(s.^2./sqrt(sum(V.^2*diag(lam),2)),1,m);
    end
    v{1}=(mean(RF,x));
    for k=1:m
        v{k+1}=(V(:,k));
        l{k+1}=getclassin('RANDVAR',varargin,RVNORMAL());
        if abs(std(l{k+1})-1)>eps
            error('les variables aleatoires doivent avoir une variance unitaire')
        end
    end
    
    
    if ischarin('pcg',varargin)
        PC = PCMODEL(RANDVARS(l(2:end)),'order',p,'pcg');
    else
        PC = PCMODEL(RANDVARS(l(2:end)),'order',p);
    end
    l = vertcat(one(PC),PC{:});
    
    L=diag([1;sqrt(lam)]);
    
    RFpc = PCRADIALMATRIX(v,size(v{1}),L*l);
    
else
    
    
    m=varargin{1};
    
    Gam = RFGAUSSIAN(0,1,RF.correl);
    RFgpc = KL(Gam,m,S,'rescale');
    RFgpc=getvalue(RFgpc);
    p = getcharin('order',varargin,1);
    
    
    H = cell(1,m);
    if ischarin('pcg',varargin)
        H(:) = {RANDPOLY(RF.marginal)};
    else
        H(:)={POLYHERMITE()};
    end
    H = RANDPOLYS(H{:});
    PC = POLYCHAOS(H,p);
    
    RFpc = decompmatrix(PC,[],[],@mytransform,RF.marginal,RFgpc,x,PC);
    
end

RFpc = FENODEFIELD(RFpc);

function r = mytransform(xi,marginal,G,x,PCnew)

xi = transfer(RANDVARS(PCnew),RANDVARS(getPC(G)),xi);
r = randomeval(G,xi);
rtemp =  cdf(RVNORMAL(),double(r));
px = eval(marginal,x);
r = zeros(size(rtemp));
for i=1:size(rtemp,2)
    r(:,i)=icdf(px,rtemp(:,i));
end

return