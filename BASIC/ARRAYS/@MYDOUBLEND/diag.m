function D = diag(varargin)
% function D = diag(V)
% function D = diag(V,k)
% Same function as DIAG but for objects of class MYDOUBLEND

if nargin>2
    error('Too many input arguments');
elseif nargin==2
    k = varargin{2};
else
    k = 0;
end
V = varargin{1};
V = double(V);

s = size2D(V);
n = s(1); % number of lines
m = s(2); % number of columns
p = n*m;
sND = sizeND(V);
pND = prod(sND);
switch min(n,m)
    case 0
        D = double.empty([0,0,sND]);
    case 1 % convert ND row or column vector to ND matrix
        dim = p+abs(k);
        D = zeros([dim,dim,sND]);
        % Construction of the list of extraction
        head = repmat(dim^2*((1:pND)-1),p,1);
        tail = repmat((dim+1)*((1:p)-1)+1,1,pND);
        if k<=0
            tail = tail-k;
        else
            tail = tail+k*dim;
        end
        list = head(:) + tail(:);
        % Extraction
        D(list) = V(:);
    otherwise % convert ND matrix to ND column vector
        if (k>-n) && (k<m) % verification of -n < k < m
            dim = min(n,m); % size of the diagonal to extract
            if k > m-dim
                dim = m-k;
            end
            if k < -n+dim
                dim = n+k;
            end
            % dim is such that -n+dim <= k <= m-dim (dim is at least 1)
            D = zeros([dim,1,sND]);
            % Construction of the list of extraction
            head = repmat(p*((1:pND)-1),dim,1);
            tail = repmat((1:n+1:dim*n),1,pND);
            if k<=0
                tail = tail-k;
            else
                tail = tail+k*n;
            end
            list = head(:) + tail(:);
            % Extraction
            D(:) = V(list);
        else
            D = double.empty([0,1,sND]);
        end
end
D = MYDOUBLEND(D);
