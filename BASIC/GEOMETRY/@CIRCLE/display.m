function display(C)
% function display(C)

disp(' ')
disp([inputname(1) ' = (' class(C) ')'])
disp(' ')
switch C.indim
    case 2
        disp(struct('cx',C.cx,'cy',C.cy,'r',C.r,'vx',C.vx,'vy',C.vy))
    case 3
        disp(struct('cx',C.cx,'cy',C.cy,'cz',C.cz,'r',C.r,'nx',C.nx,'ny',C.ny,'nz',C.nz,'vx',C.vx,'vy',C.vy))
end
