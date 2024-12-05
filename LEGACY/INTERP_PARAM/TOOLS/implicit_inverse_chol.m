function x = implicit_inverse_chol(A)

R=chol(A);

mtimes_right = @(x) R\ ( (R') \x) ;
mtimes_left  = @(x) (R\ ( (R') \(x')) )' ;

s=size(A);
x = IMPLICIT_MATRIX(s,mtimes_right,mtimes_left);

