function ke = rigi(mat,elem,xnode,xgauss,varargin)
% function ke = rigi(mat,elem,xnode,xgauss,varargin)

k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
ke = B'*k*B;

if mat.b || mat.k2 || mat.r || mat.r2 || mat.r3
    N = calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = getparam(mat,'b');
        if isa(b,'cell')
            bx = b{1};
            by = b{2};
            bx = evalparam(setparam(mat,'b',bx),'b',elem,xnode,xgauss);
            by = evalparam(setparam(mat,'b',by),'b',elem,xnode,xgauss);
            b = [bx;by];
        else
            if nargin>4 && isa(varargin{1},'POLYCHAOS')
                b = evalparampc(mat,'b',varargin{1},elem,xnode,xgauss);
            else
                b = evalparam(mat,'b',elem,xnode,xgauss);
            end
        end
        ke = ke + N'*(b'*B);
        
        if isparam(mat,'stabilize') && getparam(mat,'stabilize')>0
            bxi = norm(b);
            he = 2/sum(abs(b'*B/bxi));
            if getparam(mat,'stabilize')==1
                Pe = (bxi*he/2/k);
                tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
            elseif getparam(mat,'stabilize')==2
                tau = he/2/bxi;
            end
            ke = ke + tau*(B'*b)*(b'*B);
            % fprintf('max(Pe)=%.2d,  min(Pe)=%.2d\n',double(max(Pe)),double(min(Pe)));
        end
    end
    
    if mat.r
        r = evalparam(mat,'r',elem,xnode,xgauss);
        ke = ke + N'*r*N;
    end
    
    if mat.k2 || mat.r2 || mat.r3
        warning('la matrice de rigidite ne tient pas compte des termes non-lineaires')
    end
    % qe = varargin{1};
    % u = N*qe;
    % k2 = evalparam(mat,'k2',elem,xnode,xgauss);
    % ke = ke + B'*k2*(u.*u)*B;
    % end
end

switch getdim(elem)
    case 1
        if isparam(mat,'S')
            ke = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            ke = ke*S;
        end
    case 2
        if isparam(mat,'DIM3')
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            ke = ke*e;
        end
end
