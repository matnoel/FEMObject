r1=.04;
r2=.02;
r0=.02;
rtip=.01;
ep = 0.02;
 
lin = 0.3;
epout = .1;
lout = .3;
lmid = .5;
epsup = .5;
h1 = .15;
h2 = .3;

Point(1) =  {0 , 0 , 0 , r1} ;
Point(2) =  {lout+lmid-ep , 0 , 0 , r1} ;
Point(3) =  {lout+lmid-ep , 1-h1 , 0 , rtip} ;
Point(4) =  {lout+lmid , 1-h1, 0 , rtip} ;
Point(5) =  {lout+lmid , 0 , 0 , r0} ;
Point(6) =  {lout+lmid+lin , 0 , 0 , r0} ;
Point(7) =  {lout+lmid+lin , 1 , 0 , r1} ;
Point(8) =  {lout+ep , 1 , 0 , r1} ;
Point(9) =  {lout+ep , h2 , 0 , rtip} ;
Point(10) =  {lout , h2 , 0 , rtip} ;
Point(11) =  {lout , 1 , 0 , r2} ;
Point(12) =  {0 , 1 , 0 , r2} ;

Line(1) =  {1 , 2} ;
Line(2) =  {2 , 3} ;
Line(3) =  {3 , 4} ;
Line(4) =  {4 , 5} ;
Line(5) =  {5 , 6} ;
Line(6) =  {6 , 7} ;
Line(7) =  {7 , 8} ;
Line(8) =  {8 , 9} ;
Line(9) =  {9 , 10} ;
Line(10) =  {10 , 11} ;
Line(11) =  {11 , 12} ;
Line(12) =  {12 , 1} ;

Line Loop(101) =  {1 , 2 , 3 , 4 , 5,6,7,8,9,10,11,12} ;

Plane Surface(1) =  {101} ;


Point(13) =  {lout , 1+epout , 0 , r2} ;
Point(14) =  {0 , 1+epout , 0 , r2} ;

Line(13) =  {11 , 13} ;
Line(14) =  {13 , 14} ;
Line(15) =  {14 , 12} ;

Line Loop(102) =  {-11,13,14,15} ;

Plane Surface(2) =  {102} ;


Point(15) =  {lout , 1+epout+epsup , 0 , r1} ;
Point(16) =  {0 , 1+epout+epsup , 0 , r1} ;

Line(16) =  {13 , 15} ;
Line(17) =  {15 , 16} ;
Line(18) =  {16 , 14} ;

Line Loop(103) =  {-14,16,17,18} ;


Plane Surface(3) =  {103} ;
