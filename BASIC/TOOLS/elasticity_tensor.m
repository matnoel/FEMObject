function [ C ] = elasticity_tensor( E,nu,dim )
% function [ C ] = elasticity_tensor( E,nu,dim )
%   Return the stiffness tensor in Voigt format

if (dim~=2) && (dim~=3)
    error('elasticity_tensor.m : wrong dim')
end
C11=E*(1-nu)/((1+nu)*(1-2*nu));
C12=E*nu/((1+nu)*(1-2*nu));
Cshear=(C11-C12)/2;

if dim==2
    C=zeros(3);
    C(1,1)=C11;
    C(2,2)=C11;
    C(2,1)=C12;
    C(1,2)=C12;
    C(3,3)=Cshear;
elseif dim==3
    C=zeros(6);
    for i=1:3
        for j=1:3
            if i==j
                C(i,j)=C11;
            else
                C(i,j)=C12;
            end
        end
    end
    for i=4:6
        C(i,i)=Cshear;
    end
end

end

