function gauss = calc_gausspoints(H,n)
% gauss = calc_gausspoints(h,n)
% calcul des coordonnées et des poids des points de gauss pour des POLYFEND

M = getM(H);
if M < size(n,2)
    error('les dimensions des points de gauss ne correspondent pas aux dimensions stochastiques')
end
if size(n,2)>1
else
    n = n.*ones(1,M);
end

gausssub = cell(1,size(n,2));
for j=1:size(n,2)
    gausssub{j} = calc_gausspoints(POLYLEGENDRE(),n(1,j));
end

coord1D = cell(1,M);
w1D = cell(1,M);
for i=1:M
    I = getparam(getpoly(H,i),'I');
    j = H.j(:,i);
    nb = size(I,1);
    coord1D{i} = repmat((I(:,1)+I(:,2))/2,1,n(1,i))+(repmat(gausssub{i}.coord',nb,1).*repmat((I(:,2)-I(:,1))/2,1,n(1,i)));
    w1D{i} = repmat(gausssub{i}.w,nb,1).*(repmat(I(:,2)-I(:,1),1,n(1,i)));
    coord1D{i} = coord1D{i}(j,:);
    w1D{i} = w1D{i}(j,:);
end

nmax = max(n);
Indices = zeros(nmax^M,M);
for k=1:M
    for j=1:nmax
        for l=1:nmax^(k-1)
            Indices(l+(j-1)*(nmax)^(k-1):(nmax)^(k):end,k) = j;
        end
    end
end

if size(n,2)>1
    for k=1:M
        rep = find(Indices(:,k)>n(k));
        Indices(rep,:) = [];
    end
end
for i=1:M
    coord1D{i} = coord1D{i}(:,Indices(:,i));
    coord1D{i} = coord1D{i}';
    coord1D{i} = coord1D{i}(:);
    w1D{i} = w1D{i}(:,Indices(:,i));
    w1D{i} = w1D{i}';
    w1D{i} = w1D{i}(:);
end
gauss.coord = [coord1D{1}];
for i=2:M
    gauss.coord = [gauss.coord coord1D{i}];
end

gauss.w = ones(size(w1D{1},1),1);
for i=1:M
    gauss.w = gauss.w.*w1D{i};
end
gauss.nbgauss = length(gauss.w);
gauss.state = zeros(gauss.nbgauss,1);

gauss.state = repmat(H.state(:),1,prod(n))';
gauss.state = gauss.state(:);
