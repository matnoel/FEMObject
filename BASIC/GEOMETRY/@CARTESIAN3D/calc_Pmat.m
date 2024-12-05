function [Pe,Ps] = calc_Pmat(syscoord)
% function [Pe,Ps] = calc_Pmat(syscoord)
% Change of coordinates matrices or transition matrices Pe for strain vector epsilon
% and Ps for stress vector sigma from material coordinate system syscoord to global coordinate system
% epsilon = Pe * epsilon_glob and sigma = Ps * sigma_glob
% epsilon_glob = Ps' * epsilon and sigma_glob = Pe' * sigma

P = getbase(syscoord)'; % change of basis matrix (or transition matrix) from material coordinate system syscoord=(n1,n2,n3) to global coordinate system (x,y,z)
D1 = P.^2;
D2 = [P(2,2)*P(3,3)+P(3,2)*P(2,3) P(2,1)*P(3,3)+P(3,1)*P(2,3) P(2,1)*P(3,2)+P(3,1)*P(2,2)
    P(1,2)*P(3,3)+P(3,2)*P(1,3) P(1,1)*P(3,3)+P(3,1)*P(1,3) P(1,1)*P(3,2)+P(3,1)*P(1,2)
    P(1,2)*P(2,3)+P(2,2)*P(1,3) P(1,1)*P(2,3)+P(2,1)*P(1,3) P(1,1)*P(2,2)+P(2,1)*P(1,2)];
A = [P(2,1)*P(3,1) P(2,2)*P(3,2) P(2,3)*P(3,3)
    P(1,1)*P(3,1) P(1,2)*P(3,2) P(1,3)*P(3,3)
    P(1,1)*P(2,1) P(1,2)*P(2,2) P(1,3)*P(2,3)];
B = [P(1,2)*P(1,3) P(1,1)*P(1,3) P(1,1)*P(1,2)
    P(2,2)*P(2,3) P(2,1)*P(2,3) P(2,1)*P(2,2)
    P(3,2)*P(3,3) P(3,1)*P(3,3) P(3,1)*P(3,2)];
Pe = [D1 B; 2*A D2];
if nargout>1
    Ps = [D1 2*B; A D2];
end

end

