function B = calc_matrix2D(S,a,b,vx,vy,mx,my)
% function B = calc_matrix2D(S,a,b,vx,vy,mx,my)
% calcul de la matrice B avec
% Bij = int(li^a*lj^b)
lx = S.L;
ly = S.L;
nx = getnbpoints(lx);
ny = getnbpoints(ly);
N = nx*ny;
[I,J]=meshgrid(1:nx,1:ny);

fprintf('\ncalcul de la matrice ... ')
if a==0 && b==0
Bx = calc_matrix(lx,0,0);
By = calc_matrix(lx,0,0);
B = Bx(I(:),I(:)).*By(J(:),J(:));
elseif a==1 && b==1
Bx0 = calc_matrix(lx,0,0);
Bx1 = calc_matrix(lx,1,1);
By0 = calc_matrix(ly,0,0);
By1 = calc_matrix(ly,1,1);
B = Bx1(I(:),I(:)).*By0(J(:),J(:))+Bx0(I(:),I(:)).*By1(J(:),J(:));

elseif a==0 && b==1


if isa(vx,'double') && isa(vy,'double') 
Bx00 = calc_matrix(lx,0,0);
Bx01 = calc_matrix(lx,0,1);
By00 = calc_matrix(lx,0,0);
By01 = calc_matrix(lx,0,1);

B = vx*Bx01(I(:),I(:)).*By00(J(:),J(:))+...
    vy*Bx00(I(:),I(:)).*By01(J(:),J(:));
else
B = zeros(N,N);
if nargin>=4
mx = getnbpointsperelement(lx)+3;
my = getnbpointsperelement(ly)+3;
elseif nargin<6
mx = getnbpointsperelement(lx)+1;
my = getnbpointsperelement(ly)+1;
end

if isa(S.L,'POLYLAGRANGE') 

fprintf('calcul quadrature ... ')
gx = calc_gausslobattopoints(lx,mx);
gy = calc_gausslobattopoints(ly,my);
g = tensorize_quadrature_rule({gx,gy});

        vxval = vx(g.coord);
        vyval = vy(g.coord);
        vdu=cell(1,N);
        lval=cell(1,N);
fprintf('\nstep 1... ')
for i=1:N
    pourcentage(i,N)
    lval{i} = polyval(lx,I(i)-1,g.coord(:,1)).*polyval(ly,J(i)-1,g.coord(:,2));
    dxlval = dpolyval(lx,I(i)-1,g.coord(:,1)).*polyval(ly,J(i)-1,g.coord(:,2));
    dylval = polyval(lx,I(i)-1,g.coord(:,1)).*dpolyval(ly,J(i)-1,g.coord(:,2));
    vdu{i} =  vxval.*dxlval+vyval.*dylval;
end
fprintf('\nstep 2... ')

for i=1:N
    for j=1:N
        pourcentage((i-1)*N+j,N^2)
    B(i,j)=g.w*(lval{i}.*vdu{j});
    end
end


elseif isa(S.L,'POLYFELAGRANGE') %% a programmer
  B = zeros(N,N);  
   paramx = getparam(lx);
   paramy = getparam(ly);
    nx = paramx.n;
    ny = paramy.n;
    fprintf('Element (/%d) # ',nx*ny);
    for kx=1:nx
        for ky=1:ny
     fprintf('%d/',(kx-1)*ny+ky);
     % pourcentage((kx-1)*ny+ky,nx*ny)
    repx = paramx.reppoints{kx};
    repy = paramy.reppoints{ky};
    Nloc = length(repx)*length(repy);
    subpointsx = paramx.subpoints{kx};
    subpointsy = paramy.subpoints{ky};
    lkx = POLYLAGRANGE(subpointsx);
    lky = POLYLAGRANGE(subpointsy);
    [Iloc,Jloc]=meshgrid(repx,repy);
    
    gx = calc_gausslobattopoints(POLYLEGENDRE(),mx+1,getdomain(lkx));
    gy = calc_gausslobattopoints(POLYLEGENDRE(),my+1,getdomain(lky));
    g = tensorize_quadrature_rule({gx,gy});      
   
    
    vxval = vx(g.coord);
    vyval = vy(g.coord);
    vdu=cell(1,Nloc);
    lval=cell(1,Nloc);
        
    for i=1:length(repx)
        for j = 1:length(repy)         
            ii=(i-1)*length(repy)+j;
    lval{ii} = polyval(lkx,i-1,g.coord(:,1)).*polyval(lky,j-1,g.coord(:,2));
    dxlval = dpolyval(lkx,i-1,g.coord(:,1)).*polyval(lky,j-1,g.coord(:,2));
    dylval = polyval(lkx,i-1,g.coord(:,1)).*dpolyval(lky,j-1,g.coord(:,2));
    vdu{ii} =  vxval.*dxlval+vyval.*dylval;
        end
    end
    
     for i=1:Nloc
        for j=1:Nloc 
    ix = Iloc(i);
    iy = Jloc(i);
    ii=(ix-1)*getnbpoints(ly)+iy;
    jx = Iloc(j);
    jy = Jloc(j);
    jj=(jx-1)*getnbpoints(ly)+jy;
        
    B(ii,jj)=B(ii,jj)+...
        g.w*(lval{i}.*vdu{j});
        end
     end
    
            
        end
    end
   
end
end

end
fprintf('\n')
