function x = implicit_inverse_lu(A)

[L,U]=lu(A);

mtimes_right = @(x) U\ (L\x) ;
mtimes_left = @(x) ( L'\ ( (U')\x') )' ;

s=size(A);
x = IMPLICIT_MATRIX(s,mtimes_right,mtimes_left);

