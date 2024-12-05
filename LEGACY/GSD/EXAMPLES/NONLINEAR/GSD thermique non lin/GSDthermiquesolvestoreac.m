function l=GSDthermiquesolvestoreac(PC,tol,a1,a2,a3,b,varargin)

display_ = ischarin('display',varargin);
if display_
    fprintf('\n RESOLUTION PROB STOCHASTIQUE\n')
end

m = size(a1,1);
if length(varargin)>0 && ~isempty(varargin{1})
    l0 = varargin{1};
else
    l0 = zeros(m,PC);
end
l=l0;
opupdate=1;
modified=0;
for k=1:100
    res = b - a1*l;
    for h=1:m
        for i=1:m
            res = res - vectbase(m,h)*((l'*a3(:,:,h,i)*l)*l(i));
        end
    end

    nres = norm(res)/norm(b);
    if display_
        fprintf('iteration %d , residual = %.2d\n',k,nres);
    end

    if nres<tol
        break
    end

    if opupdate
        AT = a1;
        for i=1:m
            for j=1:m
                AT = AT + ...
                    matbase(m,m,i,j) *(l'*(a3(:,:,i,j) +2*permute(a3(:,j,:,i),[1,3,2,4]))*l ); 
            end
        end
    end


    dl = solve(AT,res);
    l = l+dl;

    errstagn = norm(l0-l)/norm(l+l0);
    if (errstagn<tol/10 || nres>1e-1) || modified==0
        opupdate=1;
    else
        opupdate=0;
    end

    l0=l;
end

if ~display_
    fprintf('(F(L) %d iterations)',k);
end

if nres>=tol
    warning('NON CONVERGENCE');
end


function I = matbase(m,n,i,j)
I = zeros(m,n);
I(i,j)=1;

return

function I = vectbase(m,i)
I = zeros(m,1);
I(i,1)=1;

return
