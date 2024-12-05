function gauss = gausstet4(order)
% function gauss = gausstet4(order)

switch order
    case {0,1}
        gauss.nbgauss = 1;
        gauss.coord = [1/4,1/4,1/4];
        gauss.w = 1/6;
        
    case {2}
        gauss.nbgauss = 4;
        gauss.coord = zeros(gauss.nbgauss,3);
        gauss.w = zeros(1,gauss.nbgauss);
        a = (5-sqrt(5))/20;
        b = (5+3*sqrt(5))/20;
        gauss.coord(1,:) = [a,a,a];
        gauss.coord(2,:) = [b,a,a];
        gauss.coord(3,:) = [a,b,a];
        gauss.coord(4,:) = [a,a,b];
        gauss.w(1:4) = 1/24;
        
    case {3}
        gauss.nbgauss = 5;
        gauss.coord = zeros(gauss.nbgauss,3);
        gauss.w = zeros(1,gauss.nbgauss);
        gauss.coord(1,:) = [1/4,1/4,1/4];
        a = 1/6;
        b = 1/2;
        gauss.coord(2,:) = [a,a,a];
        gauss.coord(3,:) = [b,a,a];
        gauss.coord(4,:) = [a,b,a];
        gauss.coord(5,:) = [a,a,b];
        gauss.w(1) = -2/15;
        gauss.w(2:5) = 3/40;
        
%     case {4}
%         gauss.nbgauss = 11;
%         gauss.coord = zeros(gauss.nbgauss,3);
%         gauss.w = zeros(1,gauss.nbgauss);
%         gauss.coord(1,:) = [1/4,1/4,1/4];
%         a = 11/14;
%         b = 1/14;
%         gauss.coord(2,:) = [a,b,b];
%         gauss.coord(3,:) = [b,a,b];
%         gauss.coord(4,:) = [b,b,a];
%         gauss.coord(5,:) = [b,b,b];
%         a = (1+sqrt(5/14))/4;
%         b = (1-sqrt(5/14))/4;
%         gauss.coord(6,:) = [a,a,b];
%         gauss.coord(7,:) = [a,b,a];
%         gauss.coord(8,:) = [a,b,b];
%         gauss.coord(9,:) = [b,a,a];
%         gauss.coord(10,:) = [b,a,b];
%         gauss.coord(11,:) = [b,b,a];
%         gauss.w(1) = -74/5625;
%         gauss.w(2:5) = 343/45000;
%         gauss.w(6:11) = 28/1125;
        
    case {4,5}
        gauss.nbgauss = 15;
        gauss.coord = zeros(gauss.nbgauss,3);
        gauss.w = zeros(1,gauss.nbgauss);
        gauss.coord(1,:) = [1/4,1/4,1/4];
        a = 1/6;
        b = 1/2;
        gauss.coord(2,:) = [a,a,a];
        gauss.coord(3,:) = [b,a,a];
        gauss.coord(4,:) = [a,b,a];
        gauss.coord(5,:) = [a,a,b];
        a = 1/8;
        b = 5/8;
        c = 3/8;
        gauss.coord(6,:) = [a,a,a];
        gauss.coord(7,:) = [b,a,a];
        gauss.coord(8,:) = [a,b,a];
        gauss.coord(9,:) = [a,a,b];
        gauss.coord(10,:) = [c,a,a];
        gauss.coord(11,:) = [a,c,a];
        gauss.coord(12,:) = [a,a,c];
        gauss.coord(13,:) = [c,c,a];
        gauss.coord(14,:) = [c,a,c];
        gauss.coord(15,:) = [a,c,c];
        gauss.w(1) = 2/45;
        gauss.w(2:5) = -27/280;
        gauss.w(6:15) = 16/315;
        
    case {6,7}
        gauss.nbgauss = 35;
        gauss.coord = zeros(gauss.nbgauss,3);
        gauss.w = zeros(1,gauss.nbgauss);
        gauss.coord(1,:) = [1/4,1/4,1/4];
        a = 1/6;
        b = 1/2;
        gauss.coord(2,:) = [a,a,a];
        gauss.coord(3,:) = [b,a,a];
        gauss.coord(4,:) = [a,b,a];
        gauss.coord(5,:) = [a,a,b];
        a = 1/8;
        b = 5/8;
        c = 3/8;
        gauss.coord(6,:) = [a,a,a];
        gauss.coord(7,:) = [b,a,a];
        gauss.coord(8,:) = [a,b,a];
        gauss.coord(9,:) = [a,a,b];
        gauss.coord(10,:) = [c,a,a];
        gauss.coord(11,:) = [a,c,a];
        gauss.coord(12,:) = [a,a,c];
        gauss.coord(13,:) = [c,c,a];
        gauss.coord(14,:) = [c,a,c];
        gauss.coord(15,:) = [a,c,c];
        a = 1/10;
        b = 7/10;
        c = 3/10;
        d = 1/2;
        gauss.coord(16,:) = [a,a,a];
        gauss.coord(17,:) = [b,a,a];
        gauss.coord(18,:) = [a,b,a];
        gauss.coord(19,:) = [a,a,b];
        gauss.coord(20,:) = [c,a,a];
        gauss.coord(21,:) = [a,c,a];
        gauss.coord(22,:) = [a,a,c];
        gauss.coord(23,:) = [d,a,a];
        gauss.coord(24,:) = [a,d,a];
        gauss.coord(25,:) = [a,a,d];
        gauss.coord(26,:) = [d,c,a];
        gauss.coord(27,:) = [d,a,c];
        gauss.coord(28,:) = [c,d,a];
        gauss.coord(29,:) = [c,a,d];
        gauss.coord(30,:) = [a,d,c];
        gauss.coord(31,:) = [a,c,d];
        gauss.coord(32,:) = [c,c,a];
        gauss.coord(33,:) = [c,a,c];
        gauss.coord(34,:) = [a,c,c];
        gauss.coord(35,:) = [c,c,c];
        gauss.w(1) = -8/945;
        gauss.w(2:5) = 243/4480;
        gauss.w(6:15) = -256/2835;
        gauss.w(16:35) = 3125/72576;
        
    otherwise
        gauss = tetrahedronquad(order);
end

