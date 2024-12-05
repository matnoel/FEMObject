function gauss=tensorize_quadrature_rule(gauss1D)

% function gauss=tensorize_quadrature_rule(gauss1D)
%
% gauss1D : cellules contenant les quadratures 1D
% gauss1D{i}.coord
% gauss1D{i}.w
% gauss1D{i}.nbgauss
%
M = length(gauss1D);
n=zeros(1,M);
for i=1:M
    n(i)=gauss1D{i}.nbgauss;
end
gauss.nbgauss = prod(n);

nmax=max(n);
Indices = zeros(nmax^M,M);
for k=1:M
    for j=1:nmax
        for l=1:nmax^(k-1)
            Indices(l+(j-1)*(nmax)^(k-1):(nmax)^(k):end,k)=j;
        end
    end
end
for k=1:M
    rep=find(Indices(:,k)>n(k));  
    Indices(rep,:)=[];
end

gauss.w=ones(1,gauss.nbgauss);
gauss.coord=zeros(gauss.nbgauss,M);

for k =1 : M
    gauss.coord(:,k)=gauss1D{k}.coord(Indices(:,k));
    gauss.w=gauss.w.*reshape(gauss1D{k}.w(Indices(:,k)),1,gauss.nbgauss);
end

