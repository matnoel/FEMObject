function [rep,P] = ispointin(D,P)
% function [rep,P] = ispointin(D,P)

tol = getfemobjectoptions('tolerancepoint');

PD = getvertices(D);

P1 = PD{1};
P2 = PD{2};
P3 = PD{3};
P4 = PD{4};

P1 = POINT(PD{1});
P2 = POINT(PD{2});
P3 = POINT(PD{3});
P4 = POINT(PD{4});

P = permute(POINT(P),[1,2,4,3]);

v1 = normalize(P2-P1);
v2 = normalize(P3-P2);
v3 = normalize(P4-P3);
v4 = normalize(P1-P4);

V1 = P-P1;
V2 = P-P2;
V3 = P-P3;
V4 = P-P4;

t1 = cross(v1,V1);
t2 = cross(v2,V2);
t3 = cross(v3,V3);
t4 = cross(v4,V4);

t = cross(v1,normalize(P3-P1));
t1 = dot(t1,t);
t2 = dot(t2,t);
t3 = dot(t3,t);
t4 = dot(t4,t);

switch getindim(D)
    case 2
        isin = t1>=-tol & t2>=-tol & t3>=-tol & t4>=-tol;
    case 3
        distortho = abs(dot(t,normalize(expand3D(V1))));
        isin = t1>=-tol & t2>=-tol & t3>=-tol & t4>=-tol & (distortho<=tol | isnan(distortho));
end

rep = find(isin);

if nargout==2
    P = P(rep);
end
