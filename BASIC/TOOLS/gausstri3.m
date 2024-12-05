function gauss = gausstri3(order)
% function gauss = gausstri3(order)

switch order
    case {0,1}
        gauss.nbgauss = 1;
        gauss.coord = [1/3,1/3];
        gauss.w = 1/2;
        
    case {2}
        gauss.nbgauss = 3;
        gauss.coord = zeros(gauss.nbgauss,2);
        gauss.w = zeros(1,gauss.nbgauss);
%         gauss.coord(1,:) = [1/2,0];
%         gauss.coord(2,:) = [1/2,1/2];
%         gauss.coord(3,:) = [0,1/2];
        gauss.coord(1,:) = [1/6,1/6];
        gauss.coord(2,:) = [2/3,1/6];
        gauss.coord(3,:) = [1/6,2/3];
        gauss.w(1) = 1/6;
        gauss.w(2) = 1/6;
        gauss.w(3) = 1/6;
        
    case {3}
        gauss.nbgauss = 4;
        gauss.coord = zeros(gauss.nbgauss,2);
        gauss.w = zeros(1,gauss.nbgauss);
%         gauss.coord(1,:) = [1,1]/3;
%         gauss.coord(2,:) = [2,11]/15;
%         gauss.coord(3,:) = [2,2]/15;
%         gauss.coord(4,:) = [11,2]/15;
        gauss.coord(1,:) = [1,1]/3;
        gauss.coord(2,:) = [1,1]/5;
        gauss.coord(3,:) = [3,1]/5;
        gauss.coord(4,:) = [1,3]/5;
        gauss.w(1) = -27/96;
        gauss.w(2) = 25/96;
        gauss.w(3) = 25/96;
        gauss.w(4) = 25/96;
        
    case {4}
        gauss.nbgauss = 6;
        gauss.coord = zeros(gauss.nbgauss,2);
        gauss.w = zeros(1,gauss.nbgauss);
        a = 0.445948490915965;
        b = 0.091576213509771;
        gauss.coord(1,:) = [a,a];
        gauss.coord(2,:) = [1-2*a,a];
        gauss.coord(3,:) = [a,1-2*a];
        gauss.coord(4,:) = [b,b];
        gauss.coord(5,:) = [1-2*b,b];
        gauss.coord(6,:) = [b,1-2*b];
        gauss.w(1:3) = 0.111690794839005;
        gauss.w(4:6) = 0.054975871827661;
        
    case {5}
        gauss.nbgauss = 7;
        gauss.coord = zeros(gauss.nbgauss,2);
        gauss.w = zeros(1,gauss.nbgauss);
        a = (6+sqrt(15))/21;
        b = 4/7-a;
        gauss.coord(1,:) = [1/3,1/3];
        gauss.coord(2,:) = [a,a];
        gauss.coord(3,:) = [1-2*a,a];
        gauss.coord(4,:) = [a,1-2*a];
        gauss.coord(5,:) = [b,b];
        gauss.coord(6,:) = [1-2*b,b];
        gauss.coord(7,:) = [b,1-2*b];
        gauss.w(1) = 9/80;
        gauss.w(2:4) = (155+sqrt(15))/2400;
        gauss.w(5:7) = 31/240-(155+sqrt(15))/2400;
        
    case {6,7}
        gauss.nbgauss = 12;
        gauss.coord = zeros(gauss.nbgauss,2);
        gauss.w = zeros(1,gauss.nbgauss);
        a = 0.063089014491502;
        b = 0.873821971016996;
        c = 0.249286745170910;
        d = 0.501426509658180;
        e = 0.310352451033785;
        f = 0.053145049844816;
        g = 0.636502499121399;
        gauss.coord(1,:) = [a,a];
        gauss.coord(2,:) = [b,a];
        gauss.coord(3,:) = [a,b];
        gauss.coord(4,:) = [c,c];
        gauss.coord(5,:) = [d,c];
        gauss.coord(6,:) = [c,d];
        gauss.coord(7,:) = [e,f];
        gauss.coord(8,:) = [f,e];
        gauss.coord(9,:) = [e,g];
        gauss.coord(10,:) = [g,e];
        gauss.coord(11,:) = [f,g];
        gauss.coord(12,:) = [g,f];
        gauss.w(1:3) = 0.025422453185103;
        gauss.w(4:6) = 0.058393137863189;
        gauss.w(7:12) = 0.041425537809187;
        
    case {8,9}
        gauss.nbgauss = 16;
        gauss.coord = zeros(gauss.nbgauss,2);
        gauss.w = zeros(1,gauss.nbgauss);
        a = 0.170569307751760;
        b = 0.658861384496480;
        c = 0.050547228317031;
        d = 0.898905543365938;
        e = 0.459292588292723;
        f = 0.081414823414554;
        g = 0.728492392955404;
        h = 0.263112829634638;
        i = 0.008394777409958;
        gauss.coord(1,:) = [1/3,1/3];
        gauss.coord(2,:) = [a,a];
        gauss.coord(3,:) = [b,a];
        gauss.coord(4,:) = [a,b];
        gauss.coord(5,:) = [c,c];
        gauss.coord(6,:) = [d,c];
        gauss.coord(7,:) = [c,d];
        gauss.coord(8,:) = [e,e];
        gauss.coord(9,:) = [f,e];
        gauss.coord(10,:) = [e,f];
        gauss.coord(11,:) = [g,h];
        gauss.coord(12,:) = [h,g];
        gauss.coord(13,:) = [g,i];
        gauss.coord(14,:) = [i,g];
        gauss.coord(15,:) = [h,i];
        gauss.coord(16,:) = [i,h];
        gauss.w(1) = 0.072157803838894;
        gauss.w(2:4) = 0.051608685267359;
        gauss.w(5:7) = 0.016229248811599;
        gauss.w(8:10) = 0.047545817133642;
        gauss.w(11:16) = 0.013615157087217;
        
    otherwise
        gauss = trianglequad(order);
        
end
