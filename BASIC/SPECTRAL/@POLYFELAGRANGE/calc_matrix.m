function B = calc_matrix(l,a,b,varargin)
% function B = calc_matrix(l,a,b)
% calcul de la matrice B avec
% Bij = int(li^a*lj^b)

param = getparam(l);
N = getnbpoints(l);
n = param.n;

if a==1
    funi = @dpolyval;
elseif a==0
    funi = @polyval;
end

if b==1
    funj = @dpolyval;
elseif b==0
    funj = @polyval;
end

B = zeros(N,N);

for k=1:n
    rep = param.reppoints{k};
    subpoints = param.subpoints{k};
    lk = POLYLAGRANGE(subpoints);
    g = calc_gausslobattopoints(POLYLEGENDRE(),param.m+1,getdomain(lk));
    
    for i=1:length(rep)
        for j=1:length(rep)
            B(rep(i),rep(j)) = B(rep(i),rep(j))+g.w*(funi(lk,i-1,g.coord).*funj(lk,j-1,g.coord));
        end
    end
    
end
