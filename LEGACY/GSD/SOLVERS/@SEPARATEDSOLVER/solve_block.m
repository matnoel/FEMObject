function [u,result] = solve_block(S,A,b)


if israndom(b)
    PC = getPC(b);
else
    PC = getPC(A);
end

u = PCTPMATRIXSUM(PC);
bu = b;

mmmax = getparam(S,'nbfoncmax');
mmmin = getparam(S,'nbfoncmin');
pfixmax=getparam(S,'pfixmax');
pfixtol=getparam(S,'pfixtol');
pfixstagn=getparam(S,'pfixstagn');
tol = getparam(S,'tol');

display_ = getparam(S,'display');
%pfixstagn=getparam(S,'pfixstagn');


if display_
    fprintf('       Tensor prod resolution (block). size = %d, \n',size(A,1))
end

normb = norm(b);
isconv=0;
errmm = zeros(1,mmmax);
l0 = cell(1,0);
for kk=1:mmmin-1
    l0{kk} = normalizephi0(normalizephi(rands(size(b,1),1,PC)));
end
n = size(b,1);

for mm=mmmin:mmmax
    
    for kk=mm
        l0{kk} = normalizephi0(normalizephi(rands(size(b,1),1,PC)));
    end
    l = l0;
    u0 = PCTPMATRIXSUM(PC);
    for kk=1:mm
        u0=u0+l{kk};
    end
    
    errpfix = zeros(1,pfixmax);
    
    
    for k=1:pfixmax
        normphi = cell(1,getnbdim(PC));
        for i=1:getnbdim(PC)
            for kk=1:mm
                l{kk}{i}=1;
            end
            mat=zeros(mm*getn(PC,i),mm*getn(PC,i));
            sec=zeros(mm*getn(PC,i),1);
            
            for kk=1:mm
                
                repkk = (kk-1)*getn(PC,i)+(1:getn(PC,i));
                a = expectnodimmtimes(i,b',l{kk});
                if all(isdetdim(a))
                    sec(repkk) = simplify(a)*mean(PC,i);
                else
                    sec(repkk)=  simplify(a);
                end
                
                for jj=1:mm
                    repjj = (jj-1)*getn(PC,i)+(1:getn(PC,i));
                    matkkjj = simplify(expectnodimmtimes(i,l{kk}',A,l{jj}));
                    
                    mat(repkk,repjj)=getmatrix(multimtimes(matkkjj',getmasseuni(PC,i)));
                end
            end
            
            mat = mat+1e-14*eye(mm*getn(PC,i));
            phi = mat\sec;
            for kk=1:mm
                repkk = (kk-1)*getn(PC,i)+(1:getn(PC,i));
                normphi{kk} = norm(phi(repkk));
                l{kk}{i} = phi(repkk)/normphi{kk};
                
            end
        end
        
        for kk=1:mm
            l{kk}{0} = 1;
        end
        mat=zeros(mm*n,mm*n);
        sec=zeros(mm*n,1);
        
        for kk=1:mm
            repkk = (kk-1)*n+(1:n);
            sec(repkk) = expect(b,l{kk});
            for jj=1:mm
                repjj = (jj-1)*n+(1:n);
                mat(repkk,repjj)= expect(l{kk}',A,l{jj});
            end
        end
        
        mat = mat+1e-14*eye(mm*n);
        phi0 = mat\sec;
        for kk=1:mm
            repkk = (kk-1)*n+(1:n);
            l{kk}{0} = phi0(repkk);
        end
        
        
        u = PCTPMATRIXSUM(PC);
        bu = b;
        for kk=1:mm
            u=u+l{kk};
            bu = bu - A*l{kk};
        end
        
        
        errpfix(k) = norm(bu)/norm(b);
        errstagn = norm(u-u0)/norm(u);
        %fprintf('iteration %3.d , nbfun %d : erreur = %.3e\n',k,mm,errpfix(k))
        fprintf('iteration %3.d , nbfun %d : stagnation = %.3e , erreur = %.3e\n',k,mm,errstagn,errpfix(k));
        %fprintf('iteration %d : erreur = %.3e\n',getm(u),norm(u-b)/norm(b));
        u0=u;l0=l;
        
        if errstagn<tol
            break
        end
        
        
        %if (k>=2 && (errpfix(k)/errpfix(k-1) > pfixstagn))
        if (k>=2 && errstagn < errpfix(k)/pfixstagn)
            %    fprintf('stagnation : ratio between error iterates = %.2e\n',errpfix(k)/errpfix(k-1));
            break
        end
    end
    
    %errmm(mm) = norm(u-b)/norm(b);
    errmm(mm) = errpfix(k);
    
    fprintf('nbfun %3.d : erreur = %.3e\n',getm(u),errmm(mm));
    if errmm(mm)<tol
        break
    end
end


result.error = errmm(1:getm(u));
