function e = conditional_expect(PC,i,fun,varargin)


type=getcharin('type',varargin,1);

if type==1
    PCi = getpcgroup(PC,i);
    RV = RANDVARS(PC);
    m = getM(PC);
    j = setdiff(1:m,i);

    if nargin==3
        Ns=10000;
    end

    if ischarin('mc',varargin) || ischarin('lhs',varargin)
        if ischarin('mc',varargin)    
            Ns = getcharin('mc',varargin);
            rs = random(RV,Ns,1);rs = [rs{:}];
        elseif ischarin('lhs',varargin)
            Ns = getcharin('lhs',varargin);
            rs = lhsrandom(RV,Ns,1);
        end
        gaussj.nbgauss = Ns;
        gaussj.coord = rs(:,j);
        gaussj.w = repmat(1/Ns,Ns,1);
    elseif ischarin('gauss',varargin) 
        gaussj = calc_gausspoints_uni(PC,getcharin('gauss',varargin));  
        gaussj = tensorize_quadrature_rule(gaussj(j));
        Ns = gaussj.nbgauss;
    end


    gaussi = calc_gausspoints(PCi,getorder(PCi)+3);

    X = zeros(Ns,m);
    X(:,j)=gaussj.coord;

    e=zeros(gaussi.nbgauss,1);
    for k=1:gaussi.nbgauss
        X(:,i) = gaussi.coord(k);
        e(k) = gaussj.w(:)'*fun(X);
    end

    ealpha = zeros(length(PCi),1);

    Hi = polyval(PCi,gaussi.coord);
    for k=1:length(PCi) 
        ealpha(k) = gaussi.w(:)'*(e(:).*Hi(:,k));
    end

    e = PCMATRIX(ealpha,[1,1],PCi);


else

    PCi = getpcgroup(PC,i);
    RV = RANDVARS(PC);
    m = getM(PC);
    j = setdiff(1:m,i);

    if nargin==3
        Ns=10000;
    end

    if ischarin('mc',varargin) || ischarin('lhs',varargin)
        if ischarin('mc',varargin)    
            Ns = getcharin('mc',varargin);
            rs = random(RV,Ns,1);rs = [rs{:}];
        elseif ischarin('lhs',varargin)
            Ns = getcharin('lhs',varargin);
            rs = lhsrandom(RV,Ns,1);
        end
        gaussj.nbgauss = Ns;
        gaussj.coord = rs;
        gaussj.w = repmat(1/Ns,Ns,1);
    elseif ischarin('gauss',varargin) 
        gaussj = calc_gausspoints_uni(PC,getcharin('gauss',varargin));  
        gaussj = tensorize_quadrature_rule(gaussj);
    end

    fX = fun(gaussj.coord);

    ealpha = zeros(length(PCi),1);

    Hi = polyval(PCi,gaussj.coord(:,i));
    for k=1:length(PCi) 
        ealpha(k) = gaussj.w(:)'*(fX(:).*Hi(:,k));
    end

    e = PCMATRIX(ealpha,[1,1],PCi);

end
