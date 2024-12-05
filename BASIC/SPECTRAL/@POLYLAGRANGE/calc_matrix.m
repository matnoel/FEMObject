function B = calc_matrix(l,a,b,varargin)
% function B = calc_matrix(l,a,b)
% calcul de la matrice B avec
% Bij = int(li^a*lj^b)

n = getnbpoints(l);
h = getparam(l,'orthopoly');
domain = getdomain(l);
g = calc_gausslobattopoints(h,n+1,domain);

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


B = zeros(n,n);
for i=1:n
    for j=1:n
        B(i,j) = g.w*(funi(l,i-1,g.coord).*funj(l,j-1,g.coord));
    end
end

