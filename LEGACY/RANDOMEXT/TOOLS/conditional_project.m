function e = conditional_project(PC,i,fun,Ns)
% function e = conditional_project(PC,i,fun,integ)

PCi = getpcgroup(PC,i);
RV = RANDVARS(PC);
m = getM(PC);
j = setdiff(1:m,i);

if isa(Ns,'double')
    rs = random(RV,Ns,1);
    rs = [rs{:}];
    gaussj.nbgauss = Ns;
    gaussj.coord = rs(:,j);
    gaussj.w = repmat(1/Ns,Ns,1);
elseif isa(Ns,'struct')
    gaussj=Ns;    
end
gaussi = calc_gausspoints(PCi,getorder(PCi)+1);

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

