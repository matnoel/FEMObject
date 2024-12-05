function N = Ntri3()

% N = inline('[1-xi(1)-xi(2),xi(1),xi(2)]','xi');
N = @(xi) [1-xi(1)-xi(2),xi(1),xi(2)];
