function B = calc_matrix(l,a,b,varargin)
% function B = calc_matrix(l,a,b)
% calcul de la matrice B avec
% Bij = int(li^a*lj^b)

p = getcharin('order',varargin);
g = calc_gausspoints(l,p+1);

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

B = zeros(p+1,p+1);
for i=1:p+1
    for j=1:p+1
        B(i,j) = g.w*(funi(l,i-1,g.coord).*funj(l,j-1,g.coord));
    end
end

