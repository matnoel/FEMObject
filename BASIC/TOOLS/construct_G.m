function G = construct_G(ma,mb)
% function G = construct_G(ma,mb)

Ma = [ma(1) ma(4) ma(5);
      ma(4) ma(2) ma(6);
      ma(5) ma(6) ma(3)];
Mb = [mb(1) mb(4) mb(5);
      mb(4) mb(2) mb(6);
      mb(5) mb(6) mb(3)];

G11 = 4*Ma.*Mb;
G12 = 2*[ma(1)*mb(4)+ma(4)*mb(1), ma(1)*mb(5)+ma(5)*mb(1), ma(4)*mb(5)+ma(5)*mb(4)
         ma(4)*mb(2)+ma(2)*mb(4), ma(4)*mb(6)+ma(6)*mb(4), ma(2)*mb(6)+ma(6)*mb(2)
         ma(5)*mb(6)+ma(6)*mb(5), ma(5)*mb(3)+ma(3)*mb(5), ma(6)*mb(3)+ma(3)*mb(6)];
G22 = [ma(1)*mb(2)+ma(2)*mb(1)+2*ma(4)*mb(4), ma(1)*mb(6)+ma(6)*mb(1)+ma(5)*mb(4)+ma(4)*mb(5), ma(4)*mb(6)+ma(6)*mb(4)+ma(5)*mb(2)+ma(2)*mb(5)
       ma(1)*mb(6)+ma(6)*mb(1)+ma(5)*mb(4)+ma(4)*mb(5), ma(1)*mb(3)+ma(3)*mb(1)+2*ma(5)*mb(5), ma(4)*mb(3)+ma(3)*mb(4)+ma(5)*mb(6)+ma(6)*mb(5)
       ma(4)*mb(6)+ma(6)*mb(4)+ma(5)*mb(2)+ma(2)*mb(5), ma(4)*mb(3)+ma(3)*mb(4)+ma(5)*mb(6)+ma(6)*mb(5), ma(2)*mb(3)+ma(3)*mb(2)+2*ma(6)*mb(6)];
G = [G11 G12
     G12' G22];
