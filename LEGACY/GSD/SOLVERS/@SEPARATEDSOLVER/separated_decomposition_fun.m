function [u,result,fs,rs] = separated_decomposition_fun(S,PC,fun,varargin)


if getparam(S,'onebyone')
    [u,result,fs,rs] = separated_decomposition_fun_onebyone(S,PC,fun,varargin{:})    ;
    return
end


rs = random(RANDVARS(PC),10000,1);
rs=[rs{:}];
u = PCTPMATRIXSUM(PC);
fs = fun(rs);
normref = sqrt(sum(fs.^2));


mmmax = getparam(S,'nbfoncmax');
mmmin = getparam(S,'nbfoncmin');
pfixmax=getparam(S,'pfixmax');
pfixstagn=getparam(S,'pfixstagn');
tol = getparam(S,'tol');


errmm = zeros(1,mmmax);
l0 = cell(1,0);
for kk=1:mmmin-1
    l0{kk} = normalizephi0(normalizephi(rands(size(b,1),1,PC)));
end

for mm=mmmin:mmmax
    
    for kk=1:mm
        l0{kk} = normalizephi0(normalizephi(rands(1,1,PC)));
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
            
            A=zeros(mm,mm);
            B=zeros(mm,getn(PC,i));
            
            for kk=1:mm
                a = double(conditional_expect(PC,i,@(x) fun(x).*randomeval(l{kk},x),varargin{:}));
                a = a(:);
                %expectnodimmtimes(i,b',l{kk});
                %if all(isdetdim(a))
                %B(kk,:) = simplify(a)*mean(PC,i);
                %else
                B(kk,:)=  simplify(a);
                %end
                for jj=1:mm
                    A(kk,jj)= simplify(expectnodimmtimes(i,l{kk}',l{jj}));
                end
            end
            
            A = A+1e-14*eye(mm);
            phi = A\B;
            for kk=1:mm
                normphi{kk} = norm(phi(kk,:));
                l{kk}{i} = phi(kk,:)'/normphi{kk};
            end
            
        end
        
        
        for kk=1:mm
            l{kk}{0} = 1;
        end
        A=zeros(mm,mm);
        B=zeros(mm,1);
        
        for kk=1:mm
            temp = conditional_expect(PC,i,@(x) fun(x).*randomeval(l{kk},x),varargin{:});
            B(kk,:) = expect(temp);
            for jj=1:mm
                A(kk,jj)= expect(l{kk},l{jj});
            end
        end
        
        A = A+1e-14*eye(mm);
        phi0 = A\B;
        for kk=1:mm
            l{kk}{0} = reshape(phi0(kk,:),1,1);
        end
        
        
        u = PCTPMATRIXSUM(PC);
        for kk=1:mm
            u=u+l{kk};
        end
        errpfix(k) = sqrt(sum((randomeval(u,rs)-fs).^2))/normref;
        errstagn = norm(u-u0)/norm(u);
        %fprintf('iteration %3.d , nbfun %d : erreur = %.3e\n',k,mm,errpfix(k))
        fprintf('iteration %3.d , nbfun %d : stagnation = %.3e , erreur = %.3e\n',k,mm,errstagn,errpfix(k));
        %fprintf('iteration %d : erreur = %.3e\n',getm(u),norm(u-b)/norm(b));
        u0=u;l0=l;
        
        
        if errpfix(k)<tol
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
