function [A,B,C,b] = advec_conv_lobatto(La,f,bc)
% function [l,A,B,C,b] = advec_conv_lobatto(La,f,bc)
% dt(u) - d2x(u) + dx(u) = f 
% La : POLYLAGRANGE ou POLYFELAGRANGE
% f function handle
% bc : boundary conditions 'dirichlet' 'periodic'
% 
% l : lagrange polynomials
% A dt(u) + B*u + C*u = b
% Aij = int(li*lj)
% Bij = int(li'*lj')
% Cij = int(li*lj')
% bi = int(li*f)

A = calc_matrix(La,0,0);
B = calc_matrix(La,1,1);
C = calc_matrix(La,0,1);

b = calc_vector(La,0,f);

n = size(A,1);
if nargin==3
    switch bc
    case 'dirichlet'
        rep = 2:n-1;
        A=A(rep,rep);
        B=B(rep,rep);
        C=C(rep,rep);

        b=b(rep);

    case 'periodic'
        A = periodic_matrix(A);
        B = periodic_matrix(B);
        C = periodic_matrix(C);
        b = periodic_vector(b);

    end

end

function A = periodic_matrix(A)
A(:,1) = A(:,1)+A(:,end);
A(1,:) = A(1,:) + A(end,:);
A(:,end)=[];
A(end,:)=[];

return

function b = periodic_vector(b)
b(1) = b(1) + b(end);
b(end)=[];
return
